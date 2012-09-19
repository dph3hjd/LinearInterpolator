//
//  LinearInterpolator.cpp
//  LinearInterpolator
//
//  Created by Hugh Dickinson on 9/18/12.
//  Copyright (c) 2012 Hugh Dickinson. All rights reserved.
//

#include "LinearInterpolator.h"
#include <iostream>


int main(int argc, char * argv[]){
    
    LinearInterpolator<2> li;
	
	li[0][0] = 0.0;
	li[0][1] = 0.5;
	li[0][2] = 1.0;
	li[1][0] = 0.5;
	li[1][1] = 0.75;
	li[1][2] = 1.0;
	li[2][0] = 1.0;
	li[2][1] = 0.25;
	li[2][2] = 0.1;
	
	li(1.5)(0.5);
		
    return 0;
}