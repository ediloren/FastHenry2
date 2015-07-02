/* This will read in the command line parameters */
/* send it argc and argv from main and a string which is a list of allowed
   option letters.  If the letter is followed by a ':' then it takes an arg */

#include "induct.h"
#include <string.h>
#ifdef SOLARIS
#include <sys/systeminfo.h>
#endif

#define GOOD 0
#define BAD 1

typedef struct _option {
  char op;       /* the option letter.  op = 'n' for '-n' */
  char *arg;     /* the argument arg = 'blah' for '-nblah' or '-n blah' */
  struct _option *next;
} Option;

Option *gather_opts();

ind_opts *Parse_Command_Line(argc, argv)
int argc;
char **argv;
{
  Option *opt_list;
  ind_opts *opts;
  int errflag;

  opts = (ind_opts *)Pmalloc(sizeof(ind_opts));
  opt_list =
    gather_opts(argc, argv, "s:m:p:o:l:f:g:a:d:r:hMk:t:b:c:e:i:D:x:S:R:v");

  /* opts left: */
  /* jnquwyzABCEFGHIJKLNOPQTUVWXYZ */

  /* get the default options */
  default_opts(opts);

  errflag = 0;
  for( ; opt_list != NULL && errflag == 0; opt_list = opt_list->next) {
    if (opt_list->op != '\0' && opt_list->arg != NULL)
      tolowercase(opt_list->arg);
    switch(opt_list->op) {
    case 's':
      if (strncmp(opt_list->arg,"ludecomp",6) == 0)
	opts->soln_technique = LUDECOMP;
      else if (strncmp(opt_list->arg,"iter",4) == 0)
	opts->soln_technique = ITERATIVE;
      /* else if (strncmp(opt_list->arg,"auto",4) == 0)
	opts->soln_technique = AUTO; */
      else {
	fprintf(stderr, "Unknown soln technique: %s\n",opt_list->arg);
	errflag++;
      }
      break;

    case 'm':
      if (strncmp(opt_list->arg,"direct",6) == 0) opts->mat_vect_prod = DIRECT;
      else if (strncmp(opt_list->arg,"multi",5) == 0) opts->mat_vect_prod = MULTIPOLE;
      /* else if (strncmp(opt_list->arg,"auto",4) == 0)
	 opts->mat_vect_prod = AUTO; */
      else {
	fprintf(stderr, "Unknown matrix vector product technique: %s\n",opt_list->arg);
	errflag++;
      }
      break;

    case 'p':
      if (read_on_off(opt_list->arg, &opts->precond) == FALSE) {
	if (strncmp(opt_list->arg, "loc", 3) == 0)
	  opts->precond = LOC;
	else if (strncmp(opt_list->arg, "sparse", 3) == 0)
	  opts->precond = SPARSE;
	else if (strncmp(opt_list->arg, "cube", 3) == 0)
	  opts->precond = SPARSE;
	else if (strncmp(opt_list->arg, "seg", 3) == 0)
	  opts->precond = SEGML;
	else if (strncmp(opt_list->arg, "diag", 3) == 0)
	  opts->precond = DIAGL;
	else if (strncmp(opt_list->arg, "posdef", 3) == 0)
	  opts->precond = POSDEF_LOC;
	else if (strncmp(opt_list->arg, "shells", 3) == 0)
	  opts->precond = SHELLS;
	else {
	  fprintf(stderr,"Unknown preconditioner setting: %s\n",opt_list->arg);
	  errflag++;
	}
      }
      break;

    case 'o':
      if (sscanf(opt_list->arg, "%d",&opts->order) != 1) {
	fprintf(stderr, "Unknown order: %s\n",opt_list->arg);
	errflag++;
      }
      break;

    case 'l':
      if (strncmp(opt_list->arg, "auto",4) ==0) opts->level = AUTO;
      else if (sscanf(opt_list->arg, "%d",&opts->level) != 1) {
	fprintf(stderr, "Unknown level: %s\n",opt_list->arg);
	errflag++;
      }
      break;

    case 'f':
      if (strcmp(opt_list->arg, "off") == 0) opts->makeFastCapFile = OFF;
      else if (strcmp(opt_list->arg, "simple") == 0)
	opts->makeFastCapFile |= SIMPLE;
      else if (strcmp(opt_list->arg, "refined") == 0)
	opts->makeFastCapFile |= REFINED;
      else if (strcmp(opt_list->arg, "both") == 0)
	opts->makeFastCapFile |= BOTH_FCAP;
      else if (strcmp(opt_list->arg, "hierarchy") == 0)
	opts->makeFastCapFile |= HIERARCHY;
      else {
	fprintf(stderr, "Unknown makeFastCapFile option: %s\n",opt_list->arg);
	errflag++;
      }
      break;

    case 'g':
      if (strcmp(opt_list->arg, "off") == 0) opts->gp_draw = OFF;
      else if (strcmp(opt_list->arg, "on") == 0)
	opts->gp_draw = THIN;
      else if (strcmp(opt_list->arg, "thin") == 0)
	opts->gp_draw = THIN;
      else if (strcmp(opt_list->arg, "thick") == 0)
	opts->gp_draw = THICK;
      else {
	fprintf(stderr, "Unknown ground plane draw setting: %s\n",opt_list->arg);
	errflag++;
      }
      break;

    case 'a':
      if (read_on_off(opt_list->arg, &opts->auto_refine) == FALSE) {
	fprintf(stderr, "Unknown auto_refine setting: %s\n",opt_list->arg);
	errflag++;
      }
      break;

    case 'i':
      if (sscanf(opt_list->arg, "%d",&opts->init_refine) != 1) {
	fprintf(stderr, "Unknown initial refinement level: %s\n",
		opt_list->arg);
	errflag++;
      }
      break;

    case 'r':
      if (sscanf(opt_list->arg, "%d",&opts->orderROM) != 1) {
	fprintf(stderr, "Unspecified order for reduced order model: %s\n",
		opt_list->arg);
	errflag++;
      }
      break;

    case 'h':
      Describe_Usage(argv[0]);   /* also exits */
      break;

    case 'M':
      opts->onlyROM = 1;
      break;

    case 'd':
      if (strcmp(opt_list->arg, "off") == 0)
	opts->dumpMats = OFF;
      else if (strcmp(opt_list->arg, "on") == 0)
	opts->dumpMats |= DUMP_ALL;
      else if (strcmp(opt_list->arg, "mrl") == 0)
	opts->dumpMats |= MRL;
      else if (strcmp(opt_list->arg, "mzmt") == 0)
	opts->dumpMats |= MZMt;
      else if (strcmp(opt_list->arg, "meshes") == 0)
	opts->dumpMats |= MESHES;
      else if (strcmp(opt_list->arg, "pre") == 0)
	opts->dumpMats |= PRE;
      else if (strcmp(opt_list->arg, "grids") == 0)
	opts->dumpMats |= GRIDS;
      else if (strcmp(opt_list->arg, "a") == 0)
	opts->dumpMats |= DUMP_A;
      else if (strcmp(opt_list->arg, "m") == 0)
	opts->dumpMats |= DUMP_M;
      else if (strcmp(opt_list->arg, "rl") == 0)
	opts->dumpMats |= DUMP_RL;
      else if (strcmp(opt_list->arg, "ls") == 0)
	opts->dumpMats |= DUMP_Ls;
      else {
	fprintf(stderr, "Unknown dumpMats option: %s\n",opt_list->arg);
	errflag++;
      }
      break;

    case 'k':
      if (strcmp(opt_list->arg, "matlab") == 0)
	opts->kind = MATLAB;
      else if (strcmp(opt_list->arg, "text") == 0)
	opts->kind = TEXT;
      else if (strcmp(opt_list->arg, "both") == 0)
	opts->kind |= BOTH_TYPES;
      else {
	fprintf(stderr, "Unknown kind of dump option: %s\n",opt_list->arg);
	errflag++;
      }
      break;

    case 't':
      if (sscanf(opt_list->arg, "%lf",&opts->tol) != 1) {
	fprintf(stderr, "Can't interpret tolerance: %s\n",opt_list->arg);
	errflag++;
      }
      break;

    case 'b':
      if (sscanf(opt_list->arg, "%lf",&opts->abs_tol) != 1) {
	fprintf(stderr, "Can't interpret absolute tolerance: %s\n",
		opt_list->arg);
	errflag++;
      }
      break;

    case 'c':
      if (sscanf(opt_list->arg, "%d",&opts->maxiters) != 1) {
	fprintf(stderr, "Can't interpret maximum number of iteration: %s\n",
		opt_list->arg);
	errflag++;
      }
      break;

    case 'e':
      if (strncmp(opt_list->arg, "auto",4) ==0) opts->limit = AUTO;
      else if (sscanf(opt_list->arg, "%d",&opts->limit) != 1) {
	fprintf(stderr, "Unknown limit: %s\n",opt_list->arg);
	errflag++;
      }
      if (opts->limit == 0) {
	fprintf(stderr, "Limit must be greater than 0\n");
	errflag++;
      }
      break;

    case 'D':
      if (read_on_off(opt_list->arg, &opts->debug) == FALSE) {
	fprintf(stderr, "Unknown debug setting: %s\n",opt_list->arg);
	errflag++;
      }
      break;

    case 'x':
      add_to_subset_of_columns(opt_list->arg, opts);
      break;

    case 'S':
       /* this points to a string of argv[] so shouldn't go away */
      opts->suffix = opt_list->arg;
      break;

    case 'R':
      if (sscanf(opt_list->arg, "%lf",&opts->shell_r0) != 1) {
	fprintf(stderr, "Can't interpret precond shell radius: %s\n",
                opt_list->arg);
	errflag++;
      }
      break;

    case 'v':
      opts->regurgitate = TRUE;
      break;

    case '\0':
      if (opts->fname == NULL) {
	if (opt_list->arg != NULL) {
	  opts->fname = opt_list->arg;
	  /*fprintf(stderr, "Filename after '-' taken as input file\n");*/
	}
	else {
	  opts->fname = "-";
	  /*fprintf(stderr, "Single '-' not understood\n");
	  errflag++;*/
	}
      }
      else {
	fprintf(stderr, "Two files specified: %s and %s. Specify only one.\n",
		opts->fname, opt_list->arg);
	errflag++;
      }
      break;

    default:
      fprintf(stderr, "Unknown option %c which has slipped past gather_opts\n",
	      opt_list->op);
      errflag++;
      break;

    }
  }

  if (errflag != 0)
    Describe_Usage(argv[0]);   /* also exits */

  fix_and_print_opts(opts);

  return opts;
}

Option *gather_opts(argc, argv, optstring)
int argc;
char **argv;
char *optstring;
{
  Option *opt_list = NULL, *opt;
  int len, count, takearg;

  count = 1;

  while(count < argc) {
    len = strlen(argv[count]);
    opt = (Option *)Pmalloc(sizeof(Option));
    if (argv[count][0] == '-') {
      if (len == 1) {
	opt->op = '\0';
	opt->arg = NULL;
      }
      else {
	opt->op = argv[count][1];
	if (is_in_optstring(opt->op, optstring, &takearg)) {
	  if (takearg == 1) {
	    if (len > 2) opt->arg = &argv[count][2];
	    else {
	      if(checkarg(count+1, argc, argv) == BAD) {
		fprintf(stderr, "for option %c\n",opt->op);
		Describe_Usage(argv[0]);
	      }
	      opt->arg = argv[count+1];
	      count++;
	    }
	  }
	  else {
	    if (len > 2) {
	      fprintf(stderr,"%s: option %c does not take an argument\n",
		      argv[0],opt->op);
	      Describe_Usage(argv[0]);
	    }
	    opt->arg = NULL;
	  }
	}
	else {
	  fprintf(stderr, "%s: Unknown option: %c\n",argv[0],opt->op);
	  Describe_Usage(argv[0]);
	}
      }
    }
    else {
      opt->op = '\0';
      opt->arg = argv[count];
    }
    opt->next = opt_list;
    opt_list = opt;
    count++;
  }

  return opt_list;
}

char *Pmalloc(size)
int size;
{
  char *blah;

  blah = (char *)malloc(size);

  if (blah == NULL) {
    fprintf(stderr, "Pmalloc: out of space trying to get %d bytes\n",size);
    exit(1);
  }

  return blah;
}

is_in_optstring(op, string, takearg)
char op, *string;
int *takearg;
{
  char *pos;

  pos = strchr(string, op);

  if (pos == NULL)
    return 0;
  else {
    if (pos[1] == ':') *takearg = 1;
    else *takearg = 0;
    return 1;
  }
}

checkarg(index, argc, argv)
int index, argc;
char **argv;
{

  if (index >= argc) {
    fprintf(stderr, "No more arguments ");
    return BAD;
  }

  if (argv[index][0] == '-') {
    fprintf(stderr,"Need an argument ");
    return BAD;
  }

  return GOOD;
}

Describe_Usage(name)
char *name;
{
  fprintf(stderr,
"Usage: %s [<input file>] [<Options>]\n",name);


  fprintf(stderr,
"FastHenry Version %s (%s)        see file default_opts.c for defaults\n",
          FHVERSION, FHDATE);
  fprintf(stderr,
" Options: (case insensitive arguments, n is an integer)\n \
  -                           = read from stdin \n \
  -s {ludecomp | iterative}   = Matrix solution method \n \
  -m {direct | multi}         = Matrix-vector product method \n \
  -p {on | off | loc | posdef = Preconditioner on or off. \n \
       | cube | seg | diag      on, sparse, cube are equivalent.\n \
       | shells                 loc and posdef are similar to release 1.5. \n \
                                seg and diag are similar to cube, but \n \
                                consume less memory and less time to factor.\n \
                                shells uses the current shell idea\n \
  -o n                        = where n = order of multipole expansions \n \
  -l {n | auto}               = where n = number of partitioning levels \n \
                                  or  auto = choose automatically \n \
  -f {off | simple | refined | both | hierarchy}\n \
                              = type of zbuf file to make. Will exit \n \
                                after completion. simple and refined refer to \n \
                                the 3D geometry. hierarchy refers to the \n \
                                nonuniform ground planes (only one allowed).\n \
  -g {off | on | thin | thick}= How to draw ground planes for -f option. \n \
                                off = draw outline only, on or thin = draw \n \
                                segments of plane as if infinitely thin, \n \
                                thick = draw completely (slower than thin).\n \
  -a {on | off}               = automatically refine the structure as if \n \
                                 doing multipole \n ");
  fprintf(stderr,
"  -i n                        = initially refine mesh to level n \n \
  -d {on | off | MRL | MZMt   = dump certain matrices and vectors to files.\n \
      | GRIDS | MESHES | pre    GRIDS and MESHES only dump in matlab format\n \
      | a | m | rl | ls}    \n \
  -k {matlab | text | both}   = kind of file to dump with -d\n \
  -t tol                      = where tol = tolerance for iteration error\n \
                                with respect to each vector element.\n \
  -b tol                      = where tol = tolerance with respect to max\n \
                                over all vector elements.\n \
  -c n                        = where n = maximum number of iterations\n \
  -e {n | auto}               = n = maximum # of filaments in a cube to be \n \
                                  considered exact. auto = (order+1)^2.\n \
  -D {on | off}               = debug information on or off. \n \
  -x portname                 = compute only the column in the admittance \n \
                                due to this port. Multiple -x are allowed\n \
  -S suffix                   = add this string to all output filenames\n \
  -r order                    - desired order for reduced order model\n \
  -M                          = compute reduced-order model and exit\n \
  -R radius                   = radius factor for shell preconditioner.\n \
                                Actual radius = radius*multipole_cube_length\n \
  -v                          = Regurgitate internal representation to stdout.\n \
                                Good for seeing what FastHenry thinks it read.\n \
");

  exit(0);
}

/*
main(argc, argv)
int argc;
char **argv;
{

  ind_opts *options;

  options = parse_command_line(argc, argv);

  printf("hello\n");

}
*/

read_on_off(str, on_off)
char *str;
int *on_off;
{
  if(strcmp(str, "on") == 0) *on_off = ON;
  else if (strcmp(str, "off") == 0) *on_off = OFF;
  else return FALSE;

  return TRUE;
}

add_to_subset_of_columns(str, opts)
char *str;
ind_opts *opts;
{
  strlist *oneport;

  tolowercase(str);
  oneport = (strlist *)Pmalloc(sizeof(strlist));

  oneport->str = (char *)Pmalloc(sizeof(char)*strlen(str)+1);
  strcpy(oneport->str, str);

  oneport->next = opts->portlist;

  opts->portlist = oneport;
}

fix_and_print_opts(opts)
ind_opts *opts;
{
  long clock;
  char hostname[BUFSIZ];

  fprintf(stdout, "Running FastHenry %s (%s)\n", FHVERSION, FHDATE);
  fprintf(stdout, "Source code fixes and compilation by FastFieldSolvers, 2012\n");

  time(&clock);
  fprintf(stdout, "  Date: %s", ctime(&clock));
#ifndef NO_GETHOSTNAME
#ifndef SOLARIS
  if(gethostname(hostname, BUFSIZ) != -1)
      fprintf(stdout, "  Host: %s\n", hostname);
  else fprintf(stdout, "  Host: ? (gethostname() failure)\n");
#else
  if (sysinfo(SI_HOSTNAME,hostname,BUFSIZ) != -1)
      fprintf(stdout, "  Host: %s\n", hostname);
  else fprintf(stdout, "  Host: ? (sysinfo() failure)\n");
#endif
#endif


  if (opts->makeFastCapFile != OFF) {
    printf("\nFastHenry visualization mode only.\n");
    printf("  Will produce files for zbuf program:\n");
    if (opts->makeFastCapFile & SIMPLE)
      printf("    zbuffile%s and zbuffile%s_shadings\n",
             opts->suffix,opts->suffix);
    if (opts->makeFastCapFile & REFINED)
      printf("    zbuffile2%s and zbuffile2%s_shadings\n",
             opts->suffix,opts->suffix);
    printf("\n");

    opts->mat_vect_prod = DIRECT;
  }
  else if (opts->soln_technique == LUDECOMP) {
    opts->mat_vect_prod = DIRECT;
    opts->precond = OFF;
    printf("Solution technique: LUDECOMP\n");
  }
  else {
    printf("Solution technique: ITERATIVE\n");
    if (opts->mat_vect_prod == DIRECT)
      printf("Matrix vector product method: DIRECT\n");
    else {
      printf("Matrix vector product method: MULTIPOLE\n");
      printf("  Order of expansion: %d\n",opts->order);
    }
    if (opts->precond != OFF)
      printf("Preconditioner: ON\n");
    else
      printf("Preconditioner: OFF\n");
    printf("Error tolerance: %lg\n",opts->tol);
  }

}
