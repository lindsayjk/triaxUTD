#pragma once
#include "common.h"

struct _Constants {
	_Constants();
	Scalar G, clight, clight2;
};

extern const _Constants Constants;
