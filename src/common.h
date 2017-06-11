#pragma once
#include <stdexcept>
#include <memory>

#define CATCH_GSL_ERRORS 1

#define throw_line_stringify_token(x) #x
#define throw_line_stringify(x) throw_line_stringify_token(x)
#define throw_line_line throw_line_stringify(__LINE__)
#define throw_line(x) throw std::runtime_error("\"" x "\" at " __FILE__  ":" throw_line_line)

typedef bool (*catch_gsl_error_callback)(int gsl_errno, void *cb_userinfo); // return false to ignore error

void init_gsl_error_handling__(void);
void begin_catch_gsl_errors__(const char* identifier, int ignore_gsl_errno = 0, catch_gsl_error_callback cb = nullptr, void *cb_userinfo = nullptr);
void end_catch_gsl_errors__(void);
#if CATCH_GSL_ERRORS
#define init_gsl_error_handling init_gsl_error_handling__
#define begin_catch_gsl_errors begin_catch_gsl_errors__
#define end_catch_gsl_errors end_catch_gsl_errors__
#else
#define init_gsl_error_handling(...)
#define begin_catch_gsl_errors(...)
#define end_catch_gsl_errors(...)
#endif

typedef double Scalar;

struct Vector2 {
	union {
		struct {
			Scalar x, y;
		};
		Scalar v[2];
	};
};

template<typename T>
class TypedArraybaseobj {
public:
	T* v;

public:
	TypedArraybaseobj(const TypedArraybaseobj&);
	TypedArraybaseobj(TypedArraybaseobj&&);
	TypedArraybaseobj& operator=(TypedArraybaseobj);
	~TypedArraybaseobj();

protected:
	const int size_in_elements;
	TypedArraybaseobj(int size_in_elements);

private:
	TypedArraybaseobj() : v(nullptr), size_in_elements(0) {}
};

// v is accessed from v[0] to v[len-1]
template<typename T>
class TypedArray1Dobj : public TypedArraybaseobj<T> {
public:
	TypedArray1Dobj(int len);

	int getLen() const { return len; };

protected:
	int len;
};

// v is accessed from v[0] to v[height*width-1] . Element (x,y) at: v[y*width+x] where 0 <= x <= width-1, 0 <= y <= height-1
template<typename T>
class TypedArray2Dobj : public TypedArraybaseobj<T> {
public:
	TypedArray2Dobj(int width,int height);

	int getWidth() const { return width; }
	int getHeight() const { return height; }

protected:
	int width, height;
};

// CAUTION: Only the following templated array types are implemented. Don't instantiate other templated types.
// Use the [Scalar|Vector2]Array[1|2]Dobj types to instantiate array objects
typedef TypedArray1Dobj<Scalar> ScalarArray1Dobj;
typedef TypedArray2Dobj<Scalar> ScalarArray2Dobj;
typedef TypedArray1Dobj<Vector2> Vector2Array1Dobj;
typedef TypedArray2Dobj<Vector2> Vector2Array2Dobj;

// Use the [Scalar|Vector2]Array[1|2]D types to create array objects and transfer ownership of array objects
/*
* Example code to create an array object:
* ScalarArray1D an_array(new ScalarArray1Dobj(100));
*
* The following transfers ownership of the array. an_array is no longer valid after this! Direct assignment is illegal.
* ScalarArray1D another_array = std::move(an_array);
*/
typedef std::unique_ptr<ScalarArray1Dobj> ScalarArray1D;
typedef std::unique_ptr<ScalarArray2Dobj> ScalarArray2D;
typedef std::unique_ptr<Vector2Array1Dobj> Vector2Array1D;
typedef std::unique_ptr<Vector2Array2Dobj> Vector2Array2D;

// Use the [Scalar|Vector2]Array[1|2]DRef types to pass references to array objects.
typedef std::unique_ptr<ScalarArray1Dobj>& ScalarArray1DRef;
typedef std::unique_ptr<ScalarArray2Dobj>& ScalarArray2DRef;
typedef std::unique_ptr<Vector2Array1Dobj>& Vector2Array1DRef;
typedef std::unique_ptr<Vector2Array2Dobj>& Vector2Array2DRef;

struct Value {
	enum ValType {
		ScalarType,
		IntType,
		UnsignedIntType,
		BoolType,
		CStringType
	};
	ValType valType;
	union {
		Scalar scalarVal;
		int intVal;
		unsigned int unsignedIntVal;
		bool boolVal;
		const char* cstrVal;
	};
};
