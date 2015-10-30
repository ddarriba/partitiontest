/*
 * FileUtilities.cpp
 *
 *  Created on: Sep 18, 2014
 *      Author: diego
 */

#include "FileUtilities.h"

/* includes for working with files/directories (Linux) */
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>
#include <cerrno>
#include <cstdlib>

int FileUtilities::isDirectory (const char* path)
{
  struct stat s_buf;

  if (stat (path, &s_buf))
  {
    return 0;
  }

  return S_ISDIR(s_buf.st_mode);
}

int FileUtilities::deleteDirectoryRecursive (const char* directory_name)
{

  if (!isDirectory (directory_name))
  {
    return ENOTDIR;
  }

  DIR *d = opendir (directory_name);
  size_t path_len = strlen (directory_name);
  int r = -1;

  if (d)
  {
    struct dirent *p;

    r = 0;

    while (!r && (p = readdir (d)))
    {
      int r2 = -1;
      char *buf;
      size_t len;

      /* Skip the names "." and ".." as we don't want to recurse on them. */
      if (!strcmp (p->d_name, ".") || !strcmp (p->d_name, ".."))
      {
        continue;
      }

      len = path_len + strlen (p->d_name) + 2;
      buf = (char *) malloc (len);

      if (buf)
      {
        struct stat statbuf;

        snprintf (buf, len, "%s/%s", directory_name, p->d_name);

        if (!stat (buf, &statbuf))
        {
          if (S_ISDIR(statbuf.st_mode))
          {
            r2 = deleteDirectoryRecursive (buf);
          }
          else
          {
            r2 = unlink (buf);
          }
        }

        free (buf);
      }

      r = r2;
    }

    closedir (d);
  }

  if (!r)
  {
    r = rmdir (directory_name);
  }

  return r;
}
