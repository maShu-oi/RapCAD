/*
 *   RapCAD - Rapid prototyping CAD IDE (www.rapcad.org)
 *   Copyright (C) 2010-2014 Giles Bathgate
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "qrotatemodule.h"
#include "context.h"
#include "numbervalue.h"
#include "vectorvalue.h"
#include "point.h"
#include "node/transformationnode.h"

QRotateModule::QRotateModule() : Module("qrotate")
{
	addParameter("quaternion");
}

Node* QRotateModule::evaluate(Context* ctx)
{
	TransformationNode* n=new TransformationNode();

	decimal w,x,y,z;
	VectorValue* vecValue=dynamic_cast<VectorValue*>(getParameterArgument(ctx,0));
	if(vecValue)
		vecValue->toQuaternion(w,x,y,z);

	decimal xx=x*x;
	decimal xy=x*y;
	decimal xz=x*z;
	decimal xw=x*w;

	decimal yy=y*y;
	decimal yz=y*z;
	decimal yw=y*w;

	decimal zz=z*z;
	decimal zw=z*w;

	decimal mx[16] = {
		1-2*(yy+zz),2*(xy-zw),2*(xz+yw),0,
		2*(xy+zw),1-2*(xx+zz),2*(yz-xw),0,
		2*(xz-yw),2*(yz+xw),1-2*(xx+yy),0,
		0,0,0,1
	};

	for(int i=0;i<16;i++)
		n->matrix[i]=mx[i];

	n->setChildren(ctx->getInputNodes());
	return n;
}
