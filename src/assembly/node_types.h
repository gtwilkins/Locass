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

#ifndef NODETYPES_H
#define NODETYPES_H

#include <vector>
#include <unordered_map>
#include <unordered_set>

class Node;

typedef std::vector<Node*> NodeList;
typedef std::vector<NodeList> NodeListList;
typedef std::unordered_set<Node*> NodeSet;
typedef std::vector<NodeSet> NodeSetList;
typedef std::unordered_map<Node*, NodeSet> NodeMap;
typedef std::unordered_map<Node*, std::pair<int32_t, int32_t> > NodeOffsetMap;
typedef std::unordered_map<Node*, int32_t> NodeIntMap;
typedef std::unordered_map<Node*, float> NodeFloatMap;
typedef std::vector< std::pair<Node*, int32_t> > NodeIntList;
typedef std::vector< std::tuple<Node*, int, int> > NodeIntIntList;
typedef std::vector< std::tuple<Node*, int, int, int> > NodeIntIntIntList;
typedef std::vector< std::pair<Node*, float> >  NodeFloatList;


#endif /* NODETYPES_H */

