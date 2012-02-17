#ifndef EXPFUNCTION_H
#define EXPFUNCTION_H

#include "function.h"
#include "context.h"

class ExpFunction : public Function
{
public:
	ExpFunction();
	Value* evaluate(Context*);
};

#endif // EXPFUNCTION_H
