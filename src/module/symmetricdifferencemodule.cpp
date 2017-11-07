/*
 *   RapCAD - Rapid prototyping CAD IDE (www.rapcad.org)
 *   Copyright (C) 2010-2017 Giles Bathgate
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

#include "symmetricdifferencemodule.h"
#include "context.h"
#include "node/symmetricdifferencenode.h"

SymmetricDifferenceModule::SymmetricDifferenceModule(Reporter* r) : Module(r,"symmetric_difference")
{
}

Node* SymmetricDifferenceModule::evaluate(Context* ctx)
{
	auto* d = new SymmetricDifferenceNode();
	d->setChildren(ctx->getInputNodes());
	return d;
}
