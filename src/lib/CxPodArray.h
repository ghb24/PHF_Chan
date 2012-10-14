#ifndef CX_PODARRAY_H
#define CX_PODARRAY_H

#include <string.h> // for memcpy
#include <stdlib.h> // for malloc/free
#include "CxDefs.h" // for assert only

namespace ct {

// swap two primitive values. Here such that we need not include <algorithm> here.
template<class FType>
void swap1(FType &A, FType &B){
   FType t = A; A = B; B = A;
}

// A reference to an existing raw pointer array. Essentially a struct ``pair of
// raw pointer and a size''. This class does not own the referenced data, it
// just provides some container-like interfaces to it.
template<class FType>
struct TArrayRef
{
   typedef FType *iterator;
   typedef FType const *const_iterator;
   typedef ::size_t size_type;

   iterator begin() { return m_pData; }
   iterator end() { return &m_pData[m_Size]; }
   const_iterator begin() const { return m_pData; }
   const_iterator end() const { return m_pData + m_Size; }

   FType &front() { assert(!empty()); return *m_pData; }
   FType &back() { assert(!empty()); return *(m_pData[m_Size-1]); }
   FType const &front() const { assert(!empty()); return *m_pData; }
   FType const &back() const { assert(!empty()); return *(m_pData[m_Size-1]); }

   FType &operator[] (size_type i) {
      assert(i < size());
      return m_pData[i];
   }
   FType const &operator[] (size_type i) const  {
      assert(i < size());
      return m_pData[i];
   }

   size_type size() const { return m_Size; }
   bool empty() const { return m_Size == 0; }

   // memset the entire array to 0.
   void clear_data() {
      if ( m_pData )
         memset(m_pData, 0, sizeof(FType)*size());
   };

   void swap(TArrayRef &other) {
      swap1(m_pData, other.m_pData);
      swap1(m_Size, other.m_Size);
   };

   TArrayRef()
      : m_pData(0), m_Size(0)
   {}

   TArrayRef(FType *pData, FType *pDataEnd)
      : m_pData(pData), m_Size(pDataEnd - pData)
   {
      assert(pDataEnd > pData);
   }

   TArrayRef(FType *pData, size_type nData)
      : m_pData(pData), m_Size(nData)
   {}
protected:
   FType
      // start of controlled array. may be 0 if no data is contained.
      *m_pData;
   size_type
      // count of actual data members.
      m_Size;
};


// A dynamic array of POD (``plain old data'') types that can be
// copied via a memcpy. Main point for this is that std::vector does not
// allow containing C-style arrays (e.g., double [3]), because C-style
// arrays are not assignable. Additionally, std::vector can be *very* slow
// when allocating large amounts of data, because the data is set to
// zero on resize. This class explicitly does not do that: It has RAII
// semantics, but non-explicitly touched data is just random.
//
// It is effectively a 'buffer-ptr + size' pair.
template<class FType>
struct TArray : public TArrayRef<FType>
{
   typedef TArrayRef<FType> FBase;
   typedef FType *iterator;
   typedef FType const *const_iterator;
   typedef ::size_t size_type;

   // WARNING: contrary to std::vector, this function *DOES NOT*
   // initialize (or touch, for that matter) the newly created data.
   // If you want that behavior, use the other resize function.
   void resize(size_type n) {
      reserve(n);
      this->m_Size = n;
   };

   void resize(size_type n, FType t) {
      size_type old_size = this->size();
      resize(n);
      if ( old_size < n ) {
         for ( size_type i = old_size; i < n; ++ i )
            this->m_pData[i] = t;
      }
   };

   void clear() {
      ::free(this->m_pData);
      this->m_pData = 0;
      this->m_Size = 0;
      m_Reserved = 0;
   };

   void resize_and_clear(size_type n) {
      resize(n);
      this->clear_data();
   }

   void push_back( FType const &t ) {
      if ( this->size() + 1 > m_Reserved ) {
         reserve(2 * this->size() + 1);
      }
      this->m_pData[this->size()] = t;
      ++this->m_Size;
   };

   void pop_back() {
      assert(!this->empty());
      this->m_Size -= 1;
   }

   void reserve(size_type n) {
      if ( m_Reserved < n ) {
         FType
            *pNewData = static_cast<FType*>(::malloc(sizeof(FType) * n));
         size_type
            nSize = this->size();
         if ( nSize != 0 )
            ::memcpy(pNewData, this->m_pData, sizeof(FType) * nSize);
         ::free(this->m_pData);
         this->m_pData = pNewData;
         this->m_Size = nSize;
         m_Reserved = n;
      }
   };


   TArray()
      : FBase(), m_Reserved(0)
   {}

   TArray(TArray const &other)
      : FBase(), m_Reserved(0)
   {
      *this = other;
   };

   // WARNING: array's content not initialized with this function!
   // (with intention!)
   explicit TArray(size_type n)
      : FBase(), m_Reserved(0)
   {
      resize(n);
   }

   TArray(size_type n, FType t)
      : FBase(), m_Reserved(0)
   {
      resize(n);
      for ( size_type i = 0; i < n; ++ i )
         this->m_pData[i] = t;
   }

   ~TArray() {
      ::free(this->m_pData);
      this->m_pData = 0;
   }

   void operator = (TArray const &other) {
      resize(other.size());
      memcpy(this->m_pData, other.m_pData, sizeof(FType) * this->size());
   };

   void swap(TArray &other) {
      swap1(this->m_pData, other.m_pData);
      swap1(this->m_Size, other.m_Size);
      swap1(this->m_Reserved, other.m_Reserved);
   };

   template<class FRandomIt>
   void assign(FRandomIt begin, FRandomIt end){
      resize(end - begin);
      for ( size_type i = 0; i < this->size(); ++ i )
         this->m_pData[i] = begin[i];
   }

   template<class FRandomIt>
   void insert(iterator pos, FRandomIt begin, FRandomIt end){
      // WARNING: not checked.
      size_type num = end - begin, ipos = pos - this->begin();
      resize(num + this->size());
      for ( size_type i = this->size(); i != ipos+num; -- i )
         this->m_pData[i-1] = this->m_pData[i-1-num];
      for ( size_type i = 0; i < num; ++ i )
         this->m_pData[ipos+i] = begin[i];
   }
private:
   size_type
      // number of entries we have space reserved for.
      m_Reserved;
};

} // namespace ct

namespace std {
   template<class FType>
   void swap( ct::TArray<FType> &A, ct::TArray<FType> &B ) {
      A.swap(B);
   }
}


#endif // CX_PODARRAY_H
