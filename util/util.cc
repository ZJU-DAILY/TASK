/*

  Copyright (C) 2011 by The Department of Computer Science and
  Technology, Tsinghua University

  Redistribution of this file is permitted under the terms of
  the BSD license.

  Author   : Dong Deng
  Created  : 2014-09-05 11:45:57 
  Modified : 2014-09-09 09:55:51
  Contact  : dd11@mails.tsinghua.edu.cn

*/


#include "util.h" 

void readData(string& filename, vector<string>& recs) {
  string str;
  ifstream input(filename, ios::in);
  while (getline(input, str)) {
    for (auto i = 0; i < str.length(); i++)
      str[i] = tolower(str[i]);
    recs.push_back(str);
  }
}
