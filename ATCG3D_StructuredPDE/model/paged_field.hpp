#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>
#if defined(__unix__) || defined(__APPLE__)
#include <sys/mman.h>
#endif

namespace atcg3d::structured_pde {

// Dense storage is the published-schema default. Sparse mode reserves a
// contiguous virtual address range backed by demand-zero anonymous pages.
// Unvisited population/work pages consume no private physical storage; element
// references and deterministic row arithmetic retain their native ordering.
template <class T> class PagedField {
  public:
    using value_type = T;
    PagedField() = default;
    ~PagedField() { release(); }
    PagedField(const PagedField &) = delete;
    PagedField &operator=(const PagedField &) = delete;
    PagedField(PagedField &&other) noexcept { swap(other); }
    PagedField &operator=(PagedField &&other) noexcept {
        release();
        dense_.clear();
        size_ = 0;
        swap(other);
        return *this;
    }
    PagedField &operator=(std::vector<T> &&values) {
        if (!sparse_) {
            dense_ = std::move(values);
            size_ = dense_.size();
        } else {
            assign(values.size(), T{});
            std::copy(values.begin(), values.end(), begin());
        }
        return *this;
    }
    void set_sparse(bool value) {
        if (size_)
            throw std::logic_error(
                "storage must be selected before allocation");
#if !defined(__unix__) && !defined(__APPLE__)
        if (value)
            throw std::runtime_error("sparse zero pages require POSIX mmap");
#endif
        sparse_ = value;
    }
    bool sparse() const noexcept { return sparse_; }
    std::size_t size() const noexcept { return size_; }
    bool empty() const noexcept { return size_ == 0; }
    T *data() noexcept { return sparse_ ? mapped_ : dense_.data(); }
    const T *data() const noexcept { return sparse_ ? mapped_ : dense_.data(); }
    T *begin() noexcept { return data(); }
    T *end() noexcept { return data() + size_; }
    const T *begin() const noexcept { return data(); }
    const T *end() const noexcept { return data() + size_; }
    T &operator[](std::size_t i) noexcept { return data()[i]; }
    const T &operator[](std::size_t i) const noexcept { return data()[i]; }
    void assign(std::size_t n, T value) {
        if (!sparse_) {
            dense_.assign(n, value);
            size_ = n;
            return;
        }
        release();
        size_ = n;
#if defined(__unix__) || defined(__APPLE__)
        if (n) {
            void *p = mmap(nullptr, n * sizeof(T), PROT_READ | PROT_WRITE,
                           MAP_PRIVATE | MAP_ANON, -1, 0);
            if (p == MAP_FAILED) {
                size_ = 0;
                throw std::bad_alloc();
            }
            mapped_ = static_cast<T *>(p);
        }
#endif
        if (value != T{})
            std::fill(begin(), end(), value);
    }
    void resize(std::size_t n) {
        if (n == size_)
            return;
        if (!sparse_) {
            dense_.resize(n);
            size_ = n;
            return;
        }
        PagedField replacement;
        replacement.set_sparse(true);
        replacement.assign(n, T{});
        std::copy_n(begin(), std::min(size_, n), replacement.begin());
        swap(replacement);
    }
    void fill(T value) {
        if (sparse_ && value == T{})
            assign(size_, value);
        else
            std::fill(begin(), end(), value);
    }
    void swap(PagedField &other) noexcept {
        dense_.swap(other.dense_);
        std::swap(mapped_, other.mapped_);
        std::swap(size_, other.size_);
        std::swap(sparse_, other.sparse_);
    }

  private:
    void release() noexcept {
#if defined(__unix__) || defined(__APPLE__)
        if (mapped_)
            munmap(mapped_, size_ * sizeof(T));
#endif
        mapped_ = nullptr;
    }
    std::vector<T> dense_;
    T *mapped_{};
    std::size_t size_{};
    bool sparse_{};
};
template <class T, class U> void fill_field(PagedField<T> &field, U value) {
    field.fill(static_cast<T>(value));
}
template <class T, class U> void fill_field(std::vector<T> &field, U value) {
    std::fill(field.begin(), field.end(), value);
}
} // namespace atcg3d::structured_pde
