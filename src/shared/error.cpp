/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "error.h"
#include <iostream>

void failure( string msg )
{
    cerr << "Error: " << msg << std::endl;
    exit( EXIT_FAILURE );
}