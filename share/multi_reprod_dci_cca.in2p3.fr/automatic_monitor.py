#!/usr/bin/env python3
"""
JUNO Batch Job Manager
Automates SLURM job submission with rate limiting and monitoring
"""

import argparse
from dataclasses import dataclass, field
import json
import time
import threading
import logging
import signal
import sys
from pathlib import Path
from queue import Queue, Empty
from typing import List, Dict, Optional, Any, Protocol
import subprocess
from datetime import datetime
from enum import Enum
from abc import ABC, abstractmethod


# ===========================================================
# Configuration & Logging
# ===========================================================

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


# ===========================================================
# Data Models
# ===========================================================

class JobStatus(Enum):
    """SLURM job status enumeration"""
    PENDING = "PENDING"
    RUNNING = "RUNNING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"
    TIMEOUT = "TIMEOUT"
    CANCELLED = "CANCELLED"
    UNKNOWN = "UNKNOWN"

JOB_STATE_MAP = {
    "PENDING": JobStatus.PENDING, "RUNNING": JobStatus.RUNNING,
    "COMPLETED": JobStatus.COMPLETED, "FAILED": JobStatus.FAILED,
    "TIMEOUT": JobStatus.TIMEOUT, "CANCELLED": JobStatus.CANCELLED,
}


@dataclass
class SLURMJob:
    """Represents a single SLURM job"""
    run: int
    jobid: int
    status: JobStatus
    created_at: datetime
    updated_at: Optional[datetime] = None
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization"""
        return {
            "run": self.run,
            "jobid": self.jobid,
            "status": self.status.value,
            "created_at": self.created_at.isoformat(),
            "updated_at": self.updated_at.isoformat() if self.updated_at else None,
            "metadata": self.metadata
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'SLURMJob':
        """Create from dictionary"""
        return cls(
            run=data["run"],
            jobid=data["jobid"],
            status=JobStatus(data["status"]),
            created_at=datetime.fromisoformat(data["created_at"]),
            updated_at=datetime.fromisoformat(data["updated_at"]) if data.get("updated_at") else None,
            metadata=data.get("metadata", {})
        )


@dataclass
class SchedulerConfig:
    """Scheduler configuration"""
    max_jobs : int = 100
    check_interval : int = 60
    skip_submitted : bool = False
    retry_failed : bool = False
    max_retries : int = 3
    enable_monitoring : bool = False
    monitoring_interval : int = 300
    show_progress : bool = True
    
    def validate(self):
        """Validate configuration values"""
        if self.max_jobs <= 0:
            raise ValueError(f"max_jobs must be positive, got {self.max_jobs}")
        if self.check_interval < 0:
            raise ValueError(f"check_interval must be non-negative, got {self.check_interval}")
        if self.max_retries < 0:
            raise ValueError(f"max_retries must be non-negative, got {self.max_retries}")


# ===========================================================
# Utility Functions
# ===========================================================

def run_cmd(cmd, shell: bool = False, timeout: float = 60.0) -> str:
    """Execute a command and return output."""
    try:
        if shell:
            if not isinstance(cmd, str):
                raise ValueError("Shell mode requires cmd to be a string")
            proc = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                check=True, 
                timeout=timeout
            )
        else:
            if not isinstance(cmd, list):
                raise ValueError("Non-shell mode requires cmd to be a list")
            proc = subprocess.run(
                cmd, 
                capture_output=True, 
                text=True, 
                check=True, 
                timeout=timeout
            )
        return proc.stdout.strip()
    except subprocess.TimeoutExpired:
        logger.error(f"Command timed out after {timeout}s: {cmd}")
        raise
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed (exit {e.returncode}): {cmd}")
        if e.stderr:
            logger.error(f"stderr: {e.stderr}")
        raise


# ===========================================================
# Base Module Class
# ===========================================================

class Module(ABC):
    """Base class for all modules with argument registration"""

    @abstractmethod
    def register_options(self, parser: argparse.ArgumentParser):
        """Register command-line arguments"""
        pass

    @abstractmethod
    def init(self, args: argparse.Namespace) -> bool:
        """Initialize with parsed arguments"""
        pass


# ===========================================================
# Throttling Strategy (Strategy Pattern)
# ===========================================================

class ThrottlingStrategy(Protocol):
    """Protocol for job submission throttling strategies"""

    @abstractmethod
    def allow_submission(self, current_jobs: int) -> bool:
        """Check if submission is allowed"""
        pass


class SimpleThrottling:
    """Simple max jobs throttling"""
    
    def __init__(self, max_jobs: int):
        self.max_jobs = max_jobs
    
    def allow_submission(self, current_jobs: int) -> bool:
        return current_jobs < self.max_jobs


class AdaptiveThrottling:
    """Adaptive throttling based on target utilization"""
    
    def __init__(self, max_jobs: int, target_utilization: float = 0.9):
        self.max_jobs = max_jobs
        self.target_utilization = target_utilization
    
    def allow_submission(self, current_jobs: int) -> bool:
        utilization = current_jobs / self.max_jobs if self.max_jobs > 0 else 0
        return utilization < self.target_utilization


# ===========================================================
# Progress Bar
# ===========================================================

class ProgressBar:
    """Thread-safe progress bar for terminal output"""
    
    def __init__(self, total: int, width: int = 50, prefix: str = "Progress"):
        self.total = total
        self.current = 0
        self.width = width
        self.prefix = prefix
        self.lock = threading.Lock()
        self.start_time = time.time()
    
    def update(self, increment: int = 1):
        """Update progress"""
        with self.lock:
            self.current += increment
            self._render()
    
    def _render(self):
        """Render progress bar"""
        if self.total == 0:
            return
        
        percent = min(self.current / self.total, 1.0)
        filled = int(self.width * percent)
        bar = '█' * filled + '░' * (self.width - filled)
        
        elapsed = time.time() - self.start_time
        if self.current > 0 and self.current < self.total:
            eta = elapsed / self.current * (self.total - self.current)
            eta_str = f"{int(eta)}s"
        else:
            eta_str = "done" if self.current >= self.total else "N/A"
        
        print(f"\r{self.prefix}: |{bar}| {self.current}/{self.total} ({percent*100:.1f}%) ETA: {eta_str}   ", 
              end='', flush=True)
        
        if self.current >= self.total:
            print()


# ===========================================================
# Run List Fetcher
# ===========================================================

class RunListFetcher(Module):
    """Fetches and filters list of runs to process"""
    
    def register_options(self, parser: argparse.ArgumentParser):
        parser.add_argument("--xrootd-url", type=str, default=None, help="XRootD server URL")
        parser.add_argument("--good-run-list", type=str, default="/eos/juno/groups/DataQuality/P25A/Physics/goodrunlist_v3.6/Physics_good_run_list.txt", help="Path to the good-run-list file")
        parser.add_argument("--lower-run", type=int, default=None, help="Process only runs newer")
        parser.add_argument("--upper-run", type=int, default=None, help="Process only runs older")

    def init(self, args: argparse.Namespace) -> bool:
        self.xrootd = args.xrootd_url
        self.good_run_list = args.good_run_list
        self.lower = args.lower_run
        self.upper = args.upper_run
        return True

    def load(self) -> List[int]:
        """Load and filter run numbers"""
        logger.info("Loading run list...")

        if self.xrootd is None:
            with open(self.good_run_list) as f:
                data = f.read().splitlines()
        else:
            cmd = f"xrdfs {self.xrootd} cat {self.good_run_list}"
            data = run_cmd(cmd, shell=True, timeout=60.0).splitlines()

        runs = []
        for line in data:
            line = line.strip()
            if not line or not line.isdigit():
                continue

            r = int(line)
            if self.lower is not None and r < self.lower:
                continue
            if self.upper is not None and r > self.upper:
                continue

            runs.append(r)

        logger.info(f"Found {len(runs)} runs (range: {min(runs)}-{max(runs)})" if runs else "No runs found")
        return runs


# ===========================================================
# SLURM Interface
# ===========================================================

class SlurmInterface(Module):
    """Interface to SLURM workload manager"""

    def register_options(self, parser : argparse.ArgumentParser):
        parser.add_argument("--user", type=str, default="traymond", help="Username used to get SLURM informations")

    def init(self, args : argparse.Namespace) -> bool:
        self.user = args.user
        return True

    def running_jobs(self) -> int:
        """Get count of running/pending jobs"""
        try:
            cmd = ["squeue", "-u", self.user, "-h", "-t", "RUNNING,PENDING"]
            out = run_cmd(cmd)
            return len(out.splitlines()) if out else 0
        except subprocess.CalledProcessError:
            logger.warning("Failed to query SLURM queue")
            return 0

    def job_state(self, jobid: int) -> JobStatus:
        """Get state of a specific job"""
        try:
            cmd = ["sacct", "-j", str(jobid), "-n", "-o", "State", "-X"]
            out = run_cmd(cmd)
            if out:
                return JOB_STATE_MAP.get(out.strip(), JobStatus.UNKNOWN)
        except subprocess.CalledProcessError:
            pass
        return JobStatus.UNKNOWN

    def show_queue(self) -> str:
        """Show user's job queue"""
        try:
            return run_cmd(["squeue", "-u", self.user])
        except subprocess.CalledProcessError as e:
            return f"Error querying queue: {e}"


# ===========================================================
# Job Registry
# ===========================================================

class JobRegistry(Module):
    """Persistent storage for job metadata"""

    def register_options(self, parser : argparse.ArgumentParser):
        parser.add_argument("--registry", type=str, default="job_registry.json", help="File used to load and save jobs informations")
        parser.add_argument("--save-interval", type=int, default=10, help="Save registry every N jobs (0 = save immediately)")

    def init(self, args : argparse.Namespace) -> bool:
        self.path = Path(args.registry)
        self.save_interval = max(0, args.save_interval)
        self.lock = threading.Lock()
        self.jobs: Dict[int, Dict[int, SLURMJob]] = {}
        self.unsaved_changes = 0
        if not self._load():
            return False
        return True

    def _load(self) -> bool:
        """Load registry from disk"""
        if not self.path.exists():
            logger.info(f"{self.path} does not exists, continue by creating a new file")
            return True

        try:
            data = json.loads(self.path.read_text())
            for run_str, jobs_dict in data.items():
                run = int(run_str)
                self.jobs[run] = {}
                for jobid_str, job_data in jobs_dict.items():
                    jobid = int(jobid_str)
                    self.jobs[run][jobid] = SLURMJob.from_dict(job_data)
            logger.info(f"Loaded {len(self.jobs)} runs from registry")
            return True
        except Exception as e:
            logger.error(f"Failed to load registry: {e}")
            backup = self.path.with_suffix('.json.backup')
            self.path.rename(backup)
            return False

    def save(self, force : bool = False) -> bool:
        """Save registry to disk"""
        with self.lock:
            if not force and (self.unsaved_changes == 0 or (self.save_interval > 0 and self.unsaved_changes < self.save_interval)):
                return False

            try:
                data = {
                    str(run): {
                        str(jid): job.to_dict() 
                        for jid, job in jobs.items()
                    } 
                    for run, jobs in self.jobs.items()
                }
                tmp = self.path.with_suffix('.tmp')
                tmp.write_text(json.dumps(data, indent=4))
                tmp.replace(self.path)
                self.unsaved_changes = 0
                return True
            except Exception as e:
                logger.error(f"Save failed: {e}")
                return False

    def run_submitted(self, run: int) -> bool:
        """Check if run has already been submitted"""
        with self.lock:
            return run in self.jobs and len(self.jobs[run]) > 0

    def add_job(self, job: SLURMJob):
        """Register a new job"""
        with self.lock:
            if job.run not in self.jobs:
                self.jobs[job.run] = {}
            self.jobs[job.run][job.jobid] = job
            self.unsaved_changes += 1
            self.save()

    def update_job_status(self, run: int, jobid: int, status: JobStatus):
        with self.lock:
            if run in self.jobs and jobid in self.jobs[run]:
                self.jobs[run][jobid].status = status
                self.jobs[run][jobid].updated_at = datetime.now()
                self.unsaved_changes += 1
                self.save()

    def get_stats(self) -> Dict[str, Any]:
        with self.lock:
            all_jobs = [job for jobs in self.jobs.values() for job in jobs.values()]
            status_counts = {s.value: sum(1 for j in all_jobs if j.status == s) 
                           for s in JobStatus}
            return {
                "total_runs": len(self.jobs), 
                "total_jobs": len(all_jobs), 
                "status_breakdown": status_counts
            }


# ===========================================================
# Job Launcher
# ===========================================================

class JobLauncher(Module):
    """Submits jobs via launcher script"""

    def register_options(self, parser : argparse.ArgumentParser):
        parser.add_argument("--launcher-script", type=str, default="job_launcher.sh", help="File used to launch jobs")
        parser.add_argument("--dry-run", action="store_true", help="Don't actually submit jobs")

    def init(self, args : argparse.Namespace) -> bool:
        self.script = Path(args.launcher_script)
        self.dry_run = args.dry_run
        if not self.dry_run and not self.script.exists():
            logger.error(f"Launcher script not found: {self.script}")
            return False
        return True

    def submit_run(self, run: int) -> List[int]:
        """Submit job for a run, return list of job IDs"""
        if self.dry_run:
            logger.info(f"[DRY-RUN] Would submit run {run}")
            return [900000 + run] # Fake job ID
        
        try:
            out = run_cmd(["bash", str(self.script), "--run-number", str(run)])
            jobids = []
            for line in out.splitlines():
                if "Submitted batch job" in line:
                    try:
                        jobids.append(int(line.split()[-1]))
                    except (ValueError, IndexError):
                        logger.warning(f"Failed to parse job ID from: {line}")

            if not jobids:
                logger.warning(f"No job IDs found in output for run {run}")

            return jobids
        except Exception as e:
            logger.error(f"Submit failed for run {run}: {e}")
            return []


# ===========================================================
# Job Monitor
# ===========================================================

class JobMonitor:
    """Background job completion monitoring"""
    
    def __init__(self, slurm: SlurmInterface, registry: JobRegistry, interval: int = 300):
        self.slurm = slurm
        self.registry = registry
        self.interval = interval
        self.running = False
        self.thread : threading.Thread = None
    
    def start(self):
        if self.running:
            return
        self.running = True
        self.thread = threading.Thread(target=self._loop, daemon=True)
        self.thread.start()
        logger.info("Job monitor started")
    
    def stop(self):
        self.running = False
        if self.thread:
            self.thread.join(timeout=5)
    
    def _loop(self):
        while self.running:
            try:
                with self.registry.lock:
                    jobs = [j for jobs in self.registry.jobs.values() for j in jobs.values()]
                active = [j for j in jobs if j.status in {JobStatus.PENDING, JobStatus.RUNNING, JobStatus.UNKNOWN}]
                
                for job in active:
                    status = self.slurm.job_state(job.jobid)
                    if status != job.status and status != JobStatus.UNKNOWN:
                        logger.info(f"Job {job.jobid} (run {job.run}): {job.status.value} → {status.value}")
                        self.registry.update_job_status(job.run, job.jobid, status)
            except Exception as e:
                logger.error(f"Monitor error: {e}")
            time.sleep(self.interval)


# ===========================================================
# Scheduler
# ===========================================================

class Scheduler(Module):
    """Manages job submission rate and queue with throttling"""

    def __init__(self, slurm : SlurmInterface, launcher : JobLauncher, registry : JobRegistry, throttling : ThrottlingStrategy = None):
        self.slurm = slurm
        self.launcher = launcher
        self.registry = registry
        self.throttling = throttling
        self.queue : Queue = Queue()
        self.running = False
        self.thread : threading.Thread = None
        self.failed_runs : List[int] = []
        self.retry_counts : Dict[int, int] = {}
        self.progress : ProgressBar = None
        self.config : SchedulerConfig = None

    def register_options(self, parser : argparse.ArgumentParser):
        parser.add_argument("--max-jobs", type=int, default=100, help="Maximum number of concurrent jobs")
        parser.add_argument("--check-interval", type=int, default=60, help="Seconds between queue checks")
        parser.add_argument("--skip-submitted", action="store_true", help="Skip runs that are already in registry")
        parser.add_argument("--retry-failed", action="store_true", help="Retry failed submissions at end of queue")
        parser.add_argument("--max-retries", type=int, default=3, help="Maximum retry attempts per run")
        parser.add_argument("--enable-monitoring", action="store_true", help="Enable monitoring")
        parser.add_argument("--monitoring-interval", type=int, default=300)
        parser.add_argument("--show-progress", action="store_true", default=True)

    def init(self, args : argparse.Namespace) -> bool:
        self.config = SchedulerConfig(
            max_jobs=args.max_jobs, check_interval=args.check_interval,
            skip_submitted=args.skip_submitted, retry_failed=args.retry_failed,
            max_retries=args.max_retries, enable_monitoring=args.enable_monitoring,
            monitoring_interval=args.monitoring_interval, show_progress=args.show_progress
        )
        self.config.validate()

        if self.throttling is None:
            self.throttling = SimpleThrottling(self.config.max_jobs)
        return True

    def enqueue_runs(self, runs: List[int]):
        """Add runs to submission queue"""
        skipped = 0
        for r in runs:
            if self.config.skip_submitted and self.registry.run_submitted(r):
                logger.debug(f"Skipping already-submitted run {r}")
                skipped += 1
                continue
            self.queue.put(r)
        
        total = len(runs) - skipped
        logger.info(f"Enqueued {total} runs" + (f" (skipped {skipped})" if skipped else ""))

        if self.config.show_progress and total > 0:
            self.progress = ProgressBar(total, prefix="Submission")

    def start(self):
        """Start scheduler thread"""
        if self.running:
            logger.warning("Scheduler already running")
            return

        self.running = True
        self.thread = threading.Thread(target=self._loop, daemon=False)
        self.thread.start()
        logger.info("Scheduler started")

    def stop(self):
        """Stop scheduler thread"""
        if not self.running:
            logger.info("Scheduler not running")
            return

        logger.info("Stopping scheduler...")
        self.running = False
        if self.thread:
            self.thread.join(timeout=10)
        self.registry.save(force=True)
        if self.failed_runs:
            logger.warning(f"Failed to submit {len(self.failed_runs)} runs: {self.failed_runs}")
        logger.info("Scheduler stopped")

    def _loop(self):
        """Main scheduler loop"""
        errors = 0
        while self.running:
            try:
                count = self.slurm.running_jobs()
                if not self.config.show_progress:
                    logger.info(f"Jobs: {count}/{self.config.max_jobs} | Queue: {self.queue.qsize()}")

                while self.throttling.allow_submission(count) and not self.queue.empty():
                    try:
                        run = self.queue.get_nowait()
                    except Empty:
                        break

                    logger.info(f"Submitting run {run}")
                    jobids = self.launcher.submit_run(run)

                    if jobids:
                        for jid in jobids:
                            self.registry.add_job(SLURMJob(run, jid, JobStatus.PENDING, datetime.now()))
                        logger.info(f"Submitted run {run}: job IDs {jobids}")
                        if self.progress:
                            self.progress.update()
                        errors = 0
                    else:
                        logger.error(f"Failed to submit run {run}")
                        if self.config.retry_failed:
                            retries = self.retry_counts.get(run, 0)
                            if retries < self.config.max_retries:
                                self.retry_counts[run] = retries + 1
                                self.queue.put(run)
                            else:
                                logger.error(f"Max retries exceeded for run {run}")
                                self.failed_runs.append(run)
                                if self.progress:
                                    self.progress.update()
                        else:
                            self.failed_runs.append(run)
                            if self.progress:
                                self.progress.update()

                    count = self.slurm.running_jobs()

                if self.queue.empty() and count == 0:
                    logger.info("Complete!" if not self.failed_runs else f"Done with {len(self.failed_runs)} failures")
                    break

                errors = 0
            except Exception as e:
                errors += 1
                logger.error(f"Scheduler error ({errors}/5): {e}")
                if errors >= 5:
                    logger.critical("Too many errors, stopping")
                    break
            
            if self.running:
                time.sleep(self.config.check_interval)


# ===========================================================
# Dependency Injection Container
# ===========================================================

@dataclass
class Container:
    """Central dependency injection container"""
    fetcher : RunListFetcher
    slurm : SlurmInterface
    registry : JobRegistry
    launcher : JobLauncher
    scheduler : Scheduler
    monitor : JobMonitor = None
    
    @classmethod
    def create(cls) -> 'Container':
        """Factory method to create container with dependencies"""
        fetcher = RunListFetcher()
        slurm = SlurmInterface()
        registry = JobRegistry()
        launcher = JobLauncher()
        scheduler = Scheduler(slurm, launcher, registry)
        monitor = JobMonitor(slurm, registry)
        
        return cls(fetcher, slurm, registry, launcher, scheduler, monitor)


# ===========================================================
# Interactive Command Interface
# ===========================================================

class CommandInterpreter:
    """Interactive command shell for monitoring"""

    def __init__(self, container : Container):
        self.container = container
        self.commands = {
            "info": (self.cmd_info, "Show status"),
            "queue": (self.cmd_queue, "Show SLURM queue"),
            "stats": (self.cmd_stats, "Show statistics"),
            "stop": (self.cmd_stop, "Stop scheduler"),
            "quit": (self.cmd_quit, "Exit"),
            "help": (self.cmd_help, "Show help"),
        }

    def run(self):
        """Run interactive command loop"""
        print("\n" + "="*60)
        print("JUNO Batch Manager - Interactive Mode")
        print("Type 'help' for available commands")
        print("="*60 + "\n")

        while True:
            try:
                line = input("> ").strip()
                if not line:
                    continue

                parts = line.split()
                cmd = parts[0]
                args = parts[1:]

                fn = self.commands.get(cmd)
                if fn:
                    fn(*args)
                else:
                    print(f"Unknown command: {cmd}")
                    print("Type 'help' for available commands")

            except EOFError:
                print("\nExiting...")
                self.cmd_quit()
            except KeyboardInterrupt:
                print("\nUse 'quit' to exit")
            except Exception as e:
                logger.error(f"Command error: {e}", exc_info=True)

    def cmd_info(self):
        """Show current status"""
        running = self.container.slurm.running_jobs()
        stats = self.container.registry.get_stats()
        failed = len(self.container.scheduler.failed_runs)

        print(f"\n{'Status':=^40}")
        print(f"Jobs running:  {running}")
        print(f"Queue size:    {self.container.scheduler.queue.qsize()}")
        print(f"Failed runs:   {failed}")
        print(f"Total runs:    {stats['total_runs']}")
        print(f"Total jobs:    {stats['total_jobs']}")
        print("="*40 + "\n")

    def cmd_queue(self):
        """Show SLURM queue"""
        print("\n" + self.container.slurm.show_queue() + "\n")

    def cmd_stats(self):
        """Show detailed statistics"""
        stats = self.container.registry.get_stats()
        print(f"\n{'Job Statistics':=^40}")
        for status, count in stats['status_breakdown'].items():
            if count > 0:
                print(f"{status:12}: {count}")
        print("="*40 + "\n")

    def cmd_stop(self):
        """Stop scheduler"""
        if not self.container.scheduler.running:
            print("Scheduler not running")
            return
        self.container.scheduler.stop()
        if self.container.monitor:
            self.container.monitor.stop()
        print("Stopped")

    def cmd_quit(self):
        """Exit program"""
        self.cmd_stop()
        sys.exit(0)

    def cmd_help(self):
        print("\nCommands:")
        for name, (_, desc) in sorted(self.commands.items()):
            print(f"  {name:8} - {desc}")
        print()


# ===========================================================
# Batch Manager - Main Orchestrator
# ===========================================================

class BatchManager:
    """Main application class"""

    def __init__(self):
        self.container = Container.create()
        self.parser = argparse.ArgumentParser(
            description="JUNO Batch Job Manager",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        # Register all module options
        for module in [self.container.fetcher, self.container.slurm, 
                      self.container.registry, self.container.launcher, 
                      self.container.scheduler]:
            module.register_options(self.parser)

        # Global options
        self.parser.add_argument("--verbose", action="store_true", help="Enable debug logging")
        self.parser.add_argument("--non-interactive", action="store_true", help="Run without interactive shell")

        # Setup signal handlers
        signal.signal(signal.SIGINT, self._signal_handler)
        signal.signal(signal.SIGTERM, self._signal_handler)

    def _signal_handler(self, signum, frame):
        """Handle shutdown signals gracefully"""
        logger.info(f"Received signal {signum}")
        self.container.scheduler.stop()
        if self.container.monitor:
            self.container.monitor.stop()
        sys.exit(0)

    def run(self) -> bool:
        """Main entry point"""
        args = self.parser.parse_args()

        # Configure logging
        if args.verbose:
            logging.getLogger().setLevel(logging.DEBUG)

        # Initialize all modules
        for module in [
            self.container.fetcher, self.container.slurm, self.container.registry, 
            self.container.launcher, self.container.scheduler
        ]:
            if not module.init(args):
                return False

        runs = self.container.fetcher.load()
        if not runs:
            logger.error("No runs to process")
            return False

        # Start services
        self.container.scheduler.enqueue_runs(runs)
        self.container.scheduler.start()

        if self.container.scheduler.config.enable_monitoring:
            self.container.monitor.start()

        # Run mode
        if args.non_interactive:
            logger.info("Non-interactive mode (Ctrl+C to stop)")
            try:
                while self.container.scheduler.running:
                    time.sleep(1)
            except KeyboardInterrupt:
                logger.info("Interrupted")
            finally:
                self.container.scheduler.stop()
                if self.container.monitor:
                    self.container.monitor.stop()
            
            if self.container.scheduler.failed_runs:
                sys.exit(1)
        else:
            CommandInterpreter(self.container).run()


# ===========================================================
# Entry Point
# ===========================================================

if __name__ == "__main__":
    if not BatchManager().run():
        sys.exit(1)
    sys.exit(0)