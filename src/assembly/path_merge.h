/*
 * Copyright (C) 2017 Glen T. Wilkins <glen.t.wilkins@gmail.com>
 * Written by Glen T. Wilkins
 * 
 * This file is part of the Locass software package <https://github.com/gtwilkins/Locass>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef PATH_MERGE_H
#define PATH_MERGE_H

#include "types.h"
#include "path_reassembly.h"

class PathMerge
{
public:
    PathMerge( Node* fork, NodeList &altPath, NodeSet &pathSet, bool drxn );
    ~PathMerge();
    
    bool merge( PathVars &pv, NodeSet &delSet, bool drxn );
    
private:
    
    Node* branch;
    vector<SeqPathMerge*> qSeqs_, tSeqs_;
};


#endif /* PATH_MERGE_H */

