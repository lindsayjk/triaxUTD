#include <gsl/gsl_errno.h>
#include <cstring>
#include <string>
#include <vector>
#include "common.h"

// Handling GSL errors

static std::vector<std::string> gsl_error_catch_identifiers_stack;

static void gsl_error_handler(const char* reason, const char* file, int line, int gsl_errno)
{
	static char tmpstr[1024];
	sprintf(tmpstr, "GSL error %d \"%s\" at (%s:%d)\nGSL error catch identifiers:", gsl_errno, reason, file, line);
	for (int n = 1, auto i = gsl_error_catch_identifiers_stack.cbegin(); i != gsl_error_catch_identifiers_stack.cend(); i++, n++) {
		sprintf(tmpstr+strlen(tmpstr), "\n[%d] %s", n, i->c_str());
	}
	throw std::runtime_error(tmpstr);
}

void init_gsl_error_handling__(void)
{
	gsl_error_catch_identifiers_stack.clear();
}

void begin_catch_gsl_errors__(const char* identifer)
{
	gsl_error_catch_identifiers_stack.push_back(std::string(identifier));
	if (gsl_error_catch_identifiers_stack.size() == 1) {
		gsl_set_error_handler(gsl_error_handler);
	}
}

void end_catch_gsl_errors__(void)
{
	gsl_error_catch_identifiers_stack.pop_back();
	if(gsl_error_catch_identifiers_stack.empty()) gsl_set_error_handler(NULL);
}

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

