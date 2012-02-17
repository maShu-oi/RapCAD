#include "expfunction.h"
#include "numbervalue.h"
#include "math.h"

ExpFunction::ExpFunction() : Function("exp")
{
}

Value* ExpFunction::evaluate(Context* ctx)
{
	NumberValue* numVal=dynamic_cast<NumberValue*>(ctx->getArgument(0,"value"));
	if(numVal) {
		double num=numVal->getNumber();

		NumberValue* result;
		result = new NumberValue(exp(num));
		return result;
	}
	return new Value();
}
