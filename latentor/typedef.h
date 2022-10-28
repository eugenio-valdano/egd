//
//  typedef.h
//  latentor
//
//  Created by Eugenio Valdano on 2/13/19.
//  Copyright © 2019 Eugenio Valdano. All rights reserved.
//

#ifndef typedef_h
#define typedef_h




#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm> //sort, unique, count_if
#include <vector>
#include <iterator>
#include <ctime>
#include <utility> // std::pair
#include <string>
#include <sstream>
#include <random>
#include <numeric>
#include <assert.h>
#include <stdexcept>
#include <functional> // bind
#include <set>
#include <forward_list> // forward list
#include <iterator> // distance
#include <thread>
#include <limits.h>
#include <zlib.h>
#include <queue>

// =========== PARSE JSON
// from https://github.com/Tencent/rapidjson
#include "/Users/eugeniovaldano/my_utilities/rapidjson-master/include/rapidjson/document.h"     // rapidjson's DOM-style API
#include "/Users/eugeniovaldano/my_utilities/rapidjson-master/include/rapidjson/filereadstream.h"



typedef unsigned int NODE;
typedef std::pair<double, NODE> NODEP;
const double INF = 1E9;


#endif /* typedef_h */
