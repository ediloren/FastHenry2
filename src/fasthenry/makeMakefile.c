/* this makes the makefile I want by brute force since the IBMs lack an
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
