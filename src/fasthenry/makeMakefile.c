/*!\page LICENSE LICENSE
 
Copyright (C) 2003 by the Board of Trustees of Massachusetts Institute of
Technology, hereafter designated as the Copyright Owners.
 
License to use, copy, modify, sell and/or distribute this software and
its documentation for any purpose is hereby granted without royalty,
subject to the following terms and conditions:
 
1.  The above copyright notice and this permission notice must
appear in all copies of the software and related documentation.
 
2.  The names of the Copyright Owners may not be used in advertising or
publicity pertaining to distribution of the software without the specific,
prior written permission of the Copyright Owners.
 
3.  THE SOFTWARE IS PROVIDED "AS-IS" AND THE COPYRIGHT OWNERS MAKE NO
REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED, BY WAY OF EXAMPLE, BUT NOT
LIMITATION.  THE COPYRIGHT OWNERS MAKE NO REPRESENTATIONS OR WARRANTIES OF
MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF THE
SOFTWARE WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS TRADEMARKS OR OTHER
RIGHTS. THE COPYRIGHT OWNERS SHALL NOT BE LIABLE FOR ANY LIABILITY OR DAMAGES
WITH RESPECT TO ANY CLAIM BY LICENSEE OR ANY THIRD PARTY ON ACCOUNT OF, OR
ARISING FROM THE LICENSE, OR ANY SUBLICENSE OR USE OF THE SOFTWARE OR ANY
SERVICE OR SUPPORT.
 
LICENSEE shall indemnify, hold harmless and defend the Copyright Owners and
their trustees, officers, employees, students and agents against any and all
claims arising out of the exercise of any rights under this Agreement,
including, without limiting the generality of the foregoing, against any
damages, losses or liabilities whatsoever with respect to death or injury to
person or damage to property arising from or out of the possession, use, or
operation of Software or Licensed Program(s) by LICENSEE or its customers.
 
*//* this makes the makefile I want by brute force since the IBMs lack an
   important feature */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#define MAXLINE 10000

char *getoneline();

main(argc, argv)
int argc;
char *argv[];
{

  char *line[2], name[1000], nameroot[1000];
  char linebeg[2][MAXLINE];
  char dir[2][50], depend[2][1000];
  int numlines = 2;
  int skip, i;

  line[0] = linebeg[0];
  line[1] = linebeg[1];
  strcpy(line[0],getoneline(stdin));
  strcpy(line[1],getoneline(stdin));
  if (line[0] == NULL || line[1] == NULL) {
    fprintf(stderr, "First line: .o files in DIR\nSecond line: .o files in MUL\n");
    exit(1);
  }

  remove_returns(line[0]);
  remove_returns(line[1]);

  fprintf(stdout, "CFLAGS = -O -DFOUR\n");
  fprintf(stdout, "DIR = .\n");
  fprintf(stdout, "MUL = $(DIR)/MattMulti\n");
  fprintf(stdout, "HEADER = $(DIR)/induct.h $(DIR)/cmplx.h $(DIR)/resusage.h\n");
  fprintf(stdout, "MULHEAD = $(MUL)/mulStruct.h $(MUL)/mulGlobal.h $(MUL)/patran.h\n");
  strcpy(dir[0],"$(DIR)");
  strcpy(dir[1],"$(MUL)");
  strcpy(depend[0], "$(HEADER) $(MULHEAD)");
  strcpy(depend[1], "$(MULHEAD)");

  fprintf(stdout, "fasthenry:\t%s %s\n",linebeg[0], linebeg[1]);
  fprintf(stdout, "\t$(CC) -o fasthenry $(CFLAGS) %s %s -lm\n",
	  linebeg[0], linebeg[1]);

  for(i = 0;  i < numlines; i++)
    while(notblankline(line[i])) {
      if (sscanf(line[i],"%s%n",name,&skip) == 1) {
	get_root(name, nameroot);
	fprintf(stdout,"%s.o:\t%s/%s.c %s\n",
		nameroot,dir[i],nameroot,depend[i]);
	fprintf(stdout,"\t$(CC) $(CFLAGS) -c %s/%s.c\n\n",dir[i],nameroot);
      }
      else {
	fprintf(stderr,"Huh? Rest of line: %s\n",line[i]);
	exit(1);
      }
      line[i] += skip;
    }

}


char *getoneline(fp)
FILE *fp;
{
  static char line[MAXLINE] = { '\0' };
  char *retchar;

  do {
    retchar = fgets(line, MAXLINE, fp); 
  } while(retchar != NULL && !notblankline(line));

  if (retchar != NULL && strlen(line) == MAXLINE - 1) 
    fprintf(stderr,"Warning: line may be too long:\n%s\n",line);

  if (retchar == NULL)
    return NULL;
  else
    return line;
}

int notblankline(string)
char *string;
{
   while( *string!='\0' && isspace(*string))
     string++;

   if (*string == '\0') return 0;
     else return 1;
}

get_root(src, dest)
char *src, *dest;
{

  while(*src != '.' && *src != '\0')
    *dest++ = *src++;

  *dest = '\0';
  if (*src != '.') {
    fprintf(stderr, "Bad .o name\n");
  }
}

remove_returns(src)
char *src;
{
  while(*src != '\0') {
    if (*src == '\n')
      *src = ' ';
    src++;
  }
}
