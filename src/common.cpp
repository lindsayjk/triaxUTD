#include <cstring>
#include "common.h"

// TypedArraybaseobj implementation

template<typename T>
TypedArraybaseobj<T>::TypedArraybaseobj(const TypedArraybaseobj<T>& other)
	: size_in_elements(other.size_in_elements), v(new T[other.size_in_elements])
{
	memcpy(v, other.v, sizeof(T)*size_in_elements);
}

template<typename T>
TypedArraybaseobj<T>::TypedArraybaseobj(TypedArraybaseobj<T>&& other)
	: size_in_elements(other.size_in_elements), v(other.v)
{
	other.v = nullptr;
}

template<typename T>
TypedArraybaseobj<T>& TypedArraybaseobj<T>::operator=(TypedArraybaseobj<T> other)
{
	std::swap(v, other.v);
	return *this;
}

template<typename T>
TypedArraybaseobj<T>::TypedArraybaseobj(int size_in_elements)
	: size_in_elements(size_in_elements), v(new T[size_in_elements])
{
}

template<typename T>
TypedArraybaseobj<T>::~TypedArraybaseobj()
{
	if(v!=nullptr) delete v;
	v = nullptr;
}

// TypedArray1Dobj implementation

template<typename T>
TypedArray1Dobj<T>::TypedArray1Dobj(int len)
	: TypedArraybaseobj<T>(len), len(len)
{
}

// TypedArray2Dobj implementation

template<typename T>
TypedArray2Dobj<T>::TypedArray2Dobj(int width, int height)
	: TypedArraybaseobj<T>(width*height), width(width), height(height)
{
}

// Explicit template instantiation for supported array types
template class TypedArraybaseobj<Scalar>;
template class TypedArraybaseobj<Vector2>;
template class TypedArray1Dobj<Scalar>;
template class TypedArray2Dobj<Scalar>;
template class TypedArray1Dobj<Vector2>;
template class TypedArray2Dobj<Vector2>;

