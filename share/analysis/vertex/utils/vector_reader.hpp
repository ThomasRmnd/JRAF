#ifndef UTILS_VECTOR_READER_HPP_
#define UTILS_VECTOR_READER_HPP_

#include <vector>

template<typename _Tp>
class vector_reader {

public:

    typedef std::vector<_Tp>                             vector_type;
    typedef typename vector_type::value_type             value_type;
    typedef typename vector_type::size_type              size_type;
    typedef typename vector_type::difference_type        difference_type;
    typedef typename vector_type::reference              reference;
    typedef typename vector_type::const_reference        const_reference;
    typedef typename vector_type::pointer                pointer;
    typedef typename vector_type::const_pointer          const_pointer;
    typedef typename vector_type::iterator               iterator;
    typedef typename vector_type::const_iterator         const_iterator;
    typedef typename vector_type::reverse_iterator       reverse_iterator;
    typedef typename vector_type::const_reverse_iterator const_reverse_iterator;

    const_iterator begin() const noexcept { return m_data->begin(); }
    const_iterator cbegin() const noexcept { return m_data->cbegin(); }

    const_iterator end() const noexcept { return m_data->end(); }
    const_iterator cend() const noexcept { return m_data->cend(); }

    const_iterator rbegin() const noexcept { return m_data->rbegin(); }
    const_iterator crbegin() const noexcept { return m_data->crbegin(); }

    const_iterator rend() const noexcept { return m_data->rend(); }
    const_iterator crend() const noexcept { return m_data->crend(); }

    vector_type* data() noexcept { return m_data; }
    const vector_type* data() const noexcept { return m_data; }

    bool empty() const noexcept { return m_data->empty(); }
    size_type size() const noexcept { return m_data->size(); }

    const_reference operator[](size_type n) const { return m_data->operator[](n); }
    const_reference at(size_type n) const { return m_data->at(n); }

private:

    vector_type* m_data = nullptr;

};

#endif // UTILS_VECTOR_READER_HPP_