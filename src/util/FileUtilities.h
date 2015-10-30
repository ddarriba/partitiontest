/*
 * FileUtilities.h
 *
 *  Created on: Sep 18, 2014
 *      Author: diego
 */

#ifndef FILEUTILITIES_H_
#define FILEUTILITIES_H_

#include <cstdio>
#include <cstring>
#include <iostream>

class FileUtilities
{
public:

  /**
   * Check the existence of a file
   *
   * @param filename The file to check
   *
   * @return true, if the file exists
   */
  static inline bool existsFile (const std::string& filename)
  {
    if (FILE *file = fopen (filename.c_str (), "r"))
    {
      fclose (file);
      return true;
    }
    else
    {
      return false;
    }
  }

  /**
   * @brief Check if a path is a directory
   *
   * @param path The path to check
   *
   * @return true, if path is directory
   */
  static int isDirectory (const char* path);

  /**
   * @brief Recursively delete a directory
   *
   * @param directory_name the directory to delete
   *
   * @return 0, if OK
   */
  static int deleteDirectoryRecursive (const char* directory_name);
};

#endif /* FILEUTILITIES_H_ */
