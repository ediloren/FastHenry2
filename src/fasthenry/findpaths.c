/* this function finds all paths from node1 to node2 */
/* this does not search through ground planes yet */

#include "induct.h"

// Enrico
#include <string.h>

/* puts all the paths in 'paths' onto the beginning of 'front' */
PATHLIST *add_to_front(paths, front)
PATHLIST *paths, *front;
{
  PATHLIST *onepath;

  if (paths == NULL)
    return front;
  else {
    onepath = paths;
    while(onepath->next != NULL)
      onepath = onepath->next;
    onepath->next = front;
    return paths;
  }
}

/* allocation for graph searching routines */
char *Gmalloc(size)
int size;
{
  char *temp;

  temp = malloc(size);

  if (temp == NULL) {
    fprintf(stderr, "Gmalloc: out of space. %d bytes needed\n", size);
  }

  return temp;
}

/* allocates space for a new path list containing one NULL path */
PATHLIST *make_new_path()
{
  PATHLIST *apath;

  apath = (PATHLIST *)Gmalloc(sizeof(PATHLIST));
  apath->path = NULL;     /* (SPATH *)Gmalloc(sizeof(SPATH)); */
/*  apath->path->next = NULL;
  apath->path->seg = seg; */
  apath->next = NULL;

  return apath;
}

/*returns the real node (not equivalenced to anything) corresponding to node */
NODES *getrealnode(node)
NODES *node;
{
  NODES *realnode;

  realnode = node->equiv;
  while(realnode != realnode->equiv)
    realnode = realnode->equiv;

/*
  if (realnode->equiv != realnode) {
    fprintf(stderr, "Hey, node->equiv chain more than one node long!\n");
    exit(1);
  }
*/

  return realnode;
}

/* returns the the node that isn't 'node' but is connected to 'tseg'. */
/* 6/7/93 - made it always return the original node (and not getrealnode()) */
NODES *getothernode(node, tseg)
NODES *node;
seg_ptr tseg;
{

  SEGMENT *seg;
  PSEUDO_SEG *pseg;
  NODES *node1, *node0, *realnode;

  if (tseg.type == NORMAL) {
    seg = (SEGMENT *)tseg.segp;
    if (seg->node[0] == node)
      return seg->node[1];
    else if (seg->node[1] == node)
      return seg->node[0];
    else {
      node0 = getrealnode(seg->node[0]);
      node1 = getrealnode(seg->node[1]);
      realnode = getrealnode(node);
      if (node0 == node1) {
	fprintf(stderr, "getothernode: node1==node0, this cause infinite loop\n");
	exit(1);
      }
      if (node0 == realnode)
	return seg->node[1];     /* return node1; */
      else if (node1 == realnode)
	return seg->node[0];     /* return node0; */
      else {
	fprintf(stderr, "Hey, node %s isn't on seg %s\n",node->name, seg->name);
	exit(1);
      }
    }
  }
  else if(tseg.type == PSEUDO) {
    pseg = (PSEUDO_SEG *)tseg.segp;
    if (pseg->node[0] == node)
      return pseg->node[1];
    else if (pseg->node[1] == node)
      return pseg->node[0];
    else {
      node0 = getrealnode(pseg->node[0]);
      node1 = getrealnode(pseg->node[1]);
      realnode = getrealnode(node);
      if (node0 == node1) {
	fprintf(stderr, "getothernode: node1==node0, this would cause an infinite loop\n");
	exit(1);
      }
      if (node0 == realnode)
	return pseg->node[1];  /* return node1; */
      else if (node1 == realnode)
	return pseg->node[0];  /* return node0; */
      else {
	fprintf(stderr, "Hey, node %s isn't on pseudo_seg\n",node->name);
	exit(1);
      }
    }
  }
  else {
    fprintf(stderr,"Huh?  Unknown type of seg: %d\n",tseg.type);
    exit(1);
  }
}

/* is this node normal */
is_normal(node)
NODES *node;
{
  if (node->type == NORMAL)
    return TRUE;
  else
    return FALSE;
}

#if 1==0
   obsolete junk
is_normal_seg(seg)
seg_ptr seg;
{
  return (seg.type == NORMAL);
}

is_pseudo_seg(seg)
seg_ptr seg;
{
  return (seg.type == PSEUDO);
}

#endif

/* is this segment in a ground plane */
is_gp(seg)
seg_ptr seg;
{
  if (seg.type == NORMAL)
    return ( ((SEGMENT *)seg.segp)->type == GPTYPE ? 1 : 0);
  else if (seg.type == PSEUDO)
    return 0;
  else {
    fprintf(stderr, "is_gp: Unknown type of seg: %d\n", seg.type);
    exit(1);
  }
}

is_gp_node(node)
NODES *node;
{
  return (node->type == GPTYPE || node->type == GPHOLE);
}

#if 1==0
   unimplemented junk
is_extern_seg(seg_p)
seg_ptr seg_p;
{
  if (seg_p.type == PSEUDO)
    return (  ((PSEUDO_SEG *)seg_p.segp)->type == EXTERNTYPE);
  else
    return (1 == 0);
}

#endif

/* is node in the list */
is_node_in_list(node, nodelist)
NODES *node;
NPATH *nodelist;
{

  NPATH *nodep;

  nodep = nodelist;
  while(nodep != NULL && getrealnode(nodep->node) != getrealnode(node))
    nodep = nodep->next;

  return (nodep == NULL ? 0 : 1);
}

/* is node (not real node) in the list */
is_orignode_in_list(node, nodelist)
NODES *node;
NPATH *nodelist;
{

  NPATH *nodep;

  nodep = nodelist;
  while(nodep != NULL && nodep->node != node)
    nodep = nodep->next;

  return (nodep == NULL ? 0 : 1);
}

/* add node to front of list */
NPATH *add_node_to_list(node, nodelist)
NODES *node;
NPATH *nodelist;
{
  NPATH *temppath;

  temppath = (NPATH *)Gmalloc(sizeof(NPATH));
  temppath->node = node;
  temppath->next = nodelist;

  return temppath;
}

/* is gp in list */
int is_gp_in_list(gp, gplist)
GROUNDPLANE *gp;
GPLIST *gplist;
{
  GPLIST *gpl;

  gpl = gplist;
  while(gpl != NULL && gpl->gp != gp)
    gpl = gpl->next;

  return (gpl == NULL ? 0 : 1);
}

/* add gp to front of gplist */
GPLIST *add_to_gplist(gp, gplist)
GROUNDPLANE *gp;
GPLIST *gplist;
{
  GPLIST *onegpl;

  onegpl = (GPLIST *)Gmalloc(sizeof(GPLIST));
  onegpl->gp = gp;
  onegpl->next = gplist;

  return onegpl;
}

/* frees each element in nodelist */
free_nodelist(nodelist)
NPATH *nodelist;
{
  NPATH *temp;

  while(nodelist != NULL) {
    temp = nodelist;
    nodelist = nodelist->next;
    free(temp);
  }
}

#if 1==0
    unimplemented obsolete junk
free_spath(path)
SPATH *path;
{
  SPATH *temp;

  while(path != NULL) {
    temp = path;
    path = path->next;
    free(temp);
  }
}
#endif

/* put path on the front of each path in pathlist */
insert_path(path, pathlist)
SPATH *path;
PATHLIST *pathlist;
{

  SPATH *pathcopy, *pathelem;

  if (pathlist != NULL) {
    while(pathlist->next != NULL) {
      pathcopy = copypath(path);
      pathelem = lastelem(pathcopy);
      pathelem->next = pathlist->path;
      pathlist->path = pathcopy;
    }
  }
}

/* allocates and fills a copy of path */
SPATH *copypath(path)
SPATH *path;
{
  SPATH *begin = NULL, *elem, *temp;

  if (path != NULL) {
    begin = (SPATH *)Gmalloc(sizeof(SPATH));
    begin->seg = path->seg;
    begin->next = NULL;
    elem = begin;
    path = path->next;
  }

  while(path != NULL) {
    temp = (SPATH *)Gmalloc(sizeof(SPATH));
    elem->next = temp;
    elem = elem->next;
    elem->seg = path->seg;
    elem->next = NULL;
    path = path->next;
  }

  return begin;
}

/* returns the last nonNULL element in the spath */
SPATH *lastelem(path)
SPATH *path;
{
  while(path->next != NULL)
    path = path->next;

  return path;
}

/* returns a pointer to the node named name.  If it is not there, it
   checks the pseudo-list and returns a real node to which the
   pseudo-name corresponds. Otherwise an error.
   Note: The node returned may be '.equiv'alenced to something else.
         Call getrealnode(node) to get the real node
*/

NODES *get_node_from_name(name, indsys)
char *name;
SYS *indsys;
{

  NODES *node;
  PSEUDO_NODE *pnode;

  node = indsys->nodes;
  while((node != NULL) && (strcmp(name, node->name) != 0))
    node = node->next;

  if (node != NULL)
    return node;
  else {
    pnode = indsys->pseudo_nodes;
    while ((pnode != NULL) && (strcmp(name, pnode->name) != 0))
      pnode = pnode->next;
    if (pnode == NULL)
      return NULL;
    else
      return pnode->node;
  }
}

equivnodes(line, indsys)
char *line;
SYS *indsys;
{
  NPATH *nlist = NULL, *nl;
  PSEUDO_NODE *pnlist = NULL, *pn;
  NODES *node, *realnode;

  int skip, i;
  char name1[80], name2[80];

  /* skip over .equiv  I think */
  if (sscanf(line,"%*s%n",&skip) != 0) {
    printf("Hey, no fair\n");
    return 1;
  }
  line += skip;

  while(notblankline(line)) {
    if (sscanf(line,"%s%n",name1,&skip) != 1) {
      printf("equivnodes: sscanf error on equiv line\n");
      return 1;
    }
    line += skip;

    node = get_node_from_name(name1, indsys);
    if (node == NULL) {
      /* node doesn't exist, add it to pseudo node list */
      pn = create_pn(name1, NULL);
      pn->next = pnlist;
      pnlist = pn;
    }
    else
      nlist = add_node_to_list(node, nlist);
  }

  if (nlist == NULL) {
    fprintf(stderr, "equivnodes: No real nodes in equiv statement\n");
    exit(1);
  }

  /* the node which all the others are to become 'equiv'ed to */
  realnode = getrealnode(nlist->node);

  /* assign the pseudo nodes to realnode */
  for(pn = pnlist; pn != NULL; pn = pn->next)
    pn->node = realnode;

  /* make the others also equivalent */
  for(nl = nlist->next; nl != NULL; nl = nl->next)
    make_equiv(nl->node, realnode);

  free_nodelist(nlist);
  append_pnlist(pnlist, indsys);

  return 0;
}

PSEUDO_NODE *create_pn(name, node)
char *name;
NODES *node;
{
  PSEUDO_NODE *pn;


  if (name[0] != 'n') {
    fprintf(stderr, "create_pn: Invalid node name %s. First letter should be 'n'\n",
	    name);
    exit(1);
  }

  pn = (PSEUDO_NODE *)MattAlloc(1, sizeof(PSEUDO_NODE));

  pn->name = (char *)MattAlloc(strlen(name) + 1, sizeof(char));

  strcpy(pn->name, name);

  pn->node = node;
  pn->next = NULL;

  return pn;
}

make_equiv(orignode, realnode)
NODES *orignode, *realnode;
{
  SEGLIST *segl;
  NODES *node;

  node = getrealnode(orignode);

  if(is_real_node(realnode) == 0 ) {
    fprintf(stderr, "make_equiv: Internal error: Can't equiv to a nonreal node\n");
    exit(1);
  }

/* can equiv before this is connected to anything, so this isn't an error
  if (realnode->connected_segs == NULL) {
    fprintf(stderr, "Internal error: make_equiv: connected_segs == NULL\n");
    exit(1);
  }
*/

  if (node == realnode) {
    fprintf(stderr,"Warning: Trying to equiv nodes that are already equiv'ed.\n");
    fprintf(stderr,"  Nodes:  %s and %s\n",orignode->name, realnode->name);
    if (is_gp_node(orignode) || is_gp_node(realnode))
      fprintf(stderr, "  Maybe ground plane is not fine enough (two nodes in same place)\n");
    return;
  }

  /* find end of connected_segs list */
  for(segl = realnode->connected_segs; segl != NULL && segl->next != NULL;
      segl = segl->next)
    ;

  /* move connected_segs on old node to realnode */
  if (segl == NULL)
    realnode->connected_segs = node->connected_segs;
  else
    segl->next = node->connected_segs;
  node->connected_segs = NULL;
  node->equiv = realnode;

}

append_pnlist(pnlist, indsys)
PSEUDO_NODE *pnlist;
SYS *indsys;
{
  PSEUDO_NODE *end;

  if (pnlist != NULL) {
    for(end = pnlist; end->next != NULL; end = end->next)
      ;

    end->next = indsys->pseudo_nodes;
    indsys->pseudo_nodes = end;
  }
}



NODES *find_next_external(node)
NODES *node;
{
    /* find next external */
  while(node != NULL && node->to_end == NULL )
      node = node->next;

  if (node != NULL && node != getrealnode(node)) {
    fprintf(stderr, "Internal err: External node not a real node. node: %s\n",
	    node->name);
    exit(1);
  }

  return node;
}

add_to_connected_segs(node, seg, pseudo_seg)
NODES *node;
SEGMENT *seg;
PSEUDO_SEG *pseudo_seg;
{
  SEGLIST *segelem;
  NODES *realnode;

  segelem = (SEGLIST *)Gmalloc(sizeof(SEGLIST));

  segelem->original = node;
  realnode = getrealnode(node);

  if (pseudo_seg == NULL) {
    segelem->seg.type = NORMAL;
    segelem->seg.segp = (void *)seg;
  }
  else if (seg == NULL) {
    segelem->seg.type = PSEUDO;
    segelem->seg.segp = (void *)pseudo_seg;
  }
  else {
    fprintf(stderr,
	    "add_to_connected_seg: Error: seg or pseudo_seg must be NULL\n");
    exit(1);
  }

  segelem->next = realnode->connected_segs;
  realnode->connected_segs = segelem;
}

remove_from_connected_segs(node, seg, pseudo_seg)
NODES *node;
SEGMENT *seg;
PSEUDO_SEG *pseudo_seg;
{
  SEGLIST *segelem, *tempsegl;
  int count;
  NODES *realnode;
  void *vptr;

  if (seg != NULL)
    vptr = seg;
  else if (pseudo_seg != NULL)
    vptr = pseudo_seg;
  else {
    fprintf(stderr, "remove_from_connected_segs: Either seg or pseudo_seg must by NULL\n");
    exit(1);
  }

  realnode = getrealnode(node);

  segelem = realnode->connected_segs;
  if (segelem == NULL) {
    fprintf(stderr, "remove_from_  : No segs connected to node %s\n", realnode->name);
    exit(1);
  }

  /* check if the first one is the seg */
  if (segelem->seg.segp == vptr && segelem->original == node) {
    realnode->connected_segs = segelem->next;
    free(segelem);
    segelem = realnode->connected_segs;
  }
  else {
    while(segelem->next != NULL &&
	  (segelem->next->seg.segp != vptr
	   || segelem->next->original != node) )
      segelem = segelem->next;

    if (segelem->next == NULL) {
      fprintf(stderr, "remove_from_  : Couldn't remove seg from node %s\n",
	      realnode->name);
      exit(1);
    }
    else {
      /* remove it from the list */
      tempsegl = segelem->next;
      segelem->next = segelem->next->next;
      free(tempsegl);
    }
  }
  /* let's check to make sure this seg isn't in the list twice */
  while(segelem != NULL &&
	(segelem->seg.segp != vptr
	 || segelem->original != node) )
    segelem = segelem->next;
  if (segelem != NULL) {
    fprintf(stderr, "remove_from_ :What?  This seg is in node %s twice\n",
	    realnode->name);
    exit(1);
  }
}


double mag(x1,y1,z1)
double x1,y1,z1;
{
  return sqrt( x1*x1 + y1*y1 + z1*z1 );
}

double magsq(x1,y1,z1)
double x1,y1,z1;
{
  return ( x1*x1 + y1*y1 + z1*z1 );
}

double dotp(x1,y1,z1,x2,y2,z2)
double x1,y1,z1,x2,y2,z2;
{
  return x1*x2 + y1*y2 + z1*z2;
}

NODES *find_nearest_gpnode(x, y, z, gp, i, j)
double x,y,z;
GROUNDPLANE *gp;
int *i, *j;
{
  static int o1 = 0, mid = 1, o2 = 2;

  double cos1, cos2, distance;
  int count1, count2;
  NODES *node;
  double *xt = gp->x, *yt = gp->y, *zt = gp->z;

  cos1 = dotp(x - xt[mid], y - yt[mid], z - zt[mid],
	      gp->ux1, gp->uy1, gp->uz1);
  cos2 = dotp(x - xt[mid], y - yt[mid], z - zt[mid],
	      gp->ux2, gp->uy2, gp->uz2);

  count1 = cos1/gp->d1 + 0.5;
  count2 = cos2/gp->d2 + 0.5;

  if (count1 > gp->seg1 || count1 < 0 || count2 > gp->seg2 || count2 < 0) {
    fprintf(stderr, "find_nearest_gpnode: point (%lg,%lg,%lg) is outside of gp %s\n",
	    x,y,z,gp->name);
    exit(1);
  }

  node = gp->pnodes[count1][count2];
  *i = count1;
  *j = count2;

  distance = mag(x - node->x, y - node->y, z - node->z);

  if ((distance - gp->unitdiag/2.0)/gp->unitdiag/2.0 > 1e-3 ) {
    fprintf(stderr, "Warning: (%lg, %lg, %lg) is more than half a cell away from plane %s\n",
	    x, y, z, gp->name);
    fprintf(stderr, "  %lg > %lg.  Is the point really on the plane?\n",
	    distance, gp->unitdiag/2.0);

    if (dotp(gp->ux1,gp->uy1,gp->uz1,gp->ux2,gp->uy2,gp->uz2) > 1e-13) {
      fprintf(stderr,"find_nearest_gpnode: Edges of plane not perpendicular.\n");
      fprintf(stderr," This is not yet supported in this function \n");
    }
  }

  return node;
}

is_real_node(node)
NODES *node;
{
  return (node == node->equiv ? 1 : 0);
}

SPATH *add_seg_to_list(seg, seglist)
seg_ptr seg;
SPATH *seglist;
{
  SPATH *temppath;

  temppath = (SPATH *)Gmalloc(sizeof(SPATH));
  temppath->seg = seg;
  temppath->next = seglist;

  return temppath;
}

PSEUDO_SEG *make_pseudo_seg(node1, node2, type)
NODES *node1, *node2;
char type;
{
  PSEUDO_SEG *temp_seg;

  temp_seg = (PSEUDO_SEG *)Gmalloc(sizeof(PSEUDO_SEG));

  temp_seg->node[0] = node1;
  temp_seg->node[1] = node2;
  temp_seg->type = type;
  temp_seg->loops = NULL;
  temp_seg->is_deleted = 0;
  add_to_connected_segs(node1, NULL, temp_seg);
  add_to_connected_segs(node2, NULL, temp_seg);

#if 1==0
  unimplemented obsolete junk
  /* guess the upper bound for number of real segs in this pseudo seg */
  if (type == GPTYPE) {
    if (node1->gp != node2->gp) {
      fprintf(stderr,"Internal Err: gp pseudo_seg nodes not in same plane!\n");
      exit(1);
    }
    /* estimate path length as distance from opposite corners */
    temp_seg->upper_num_segs = node1->gp->seg1 + node1->gp->seg2;
  }
  else if (type == EXTERNTYPE)
    temp_seg->upper_num_segs = 0;
  else {
    fprintf(stderr, "Internal Error: Unknown type of pseudo_seg\n");
    exit(1);
  }
#endif

  return temp_seg;
}

SPATH *make_new_fake_segs(node, nodelist, seg_list)
NODES *node;
NPATH *nodelist;
SPATH *seg_list;
{
  NODES *anode;
  PSEUDO_SEG *pseg;
  SPATH *pathptr;
  seg_ptr seg;

  if (nodelist != NULL) {
    /* make a new fake seg between the node and the nearest one in the list. */
    /* This should reduce the length of gp meshes hopefully.  3/96           */
    anode = find_nearest_node(node, nodelist);
    /* anode = nodelist->node;  just use the first in the list (FH 2.0) */
    pseg = make_pseudo_seg(node, anode, GPTYPE);
    seg.segp = (void *)pseg;
    seg.type = PSEUDO;
    seg_list = add_seg_to_list(seg, seg_list);
  }

  return seg_list;
}

EXTERNAL *add_to_external_list(ex, ex_list)
EXTERNAL *ex, *ex_list;
{
  ex->next = ex_list;
  return ex;
}

EXTERNAL *make_external(source, Yindex, name1, name2, portname)
PSEUDO_SEG *source;
int Yindex;
char *name1, *name2, *portname;
{
 EXTERNAL *temp_ex;

 temp_ex = (EXTERNAL *)Gmalloc(sizeof(EXTERNAL));
 temp_ex->source = source;
 temp_ex->indices = NULL;
 temp_ex->Yindex = Yindex;
 temp_ex->loops = NULL;
 temp_ex->next = NULL;

 temp_ex->name1 = (char *)MattAlloc(strlen(name1)+1, sizeof(char));
 temp_ex->name2 = (char *)MattAlloc(strlen(name2)+1, sizeof(char));
 temp_ex->portname = (char *)MattAlloc(strlen(portname)+1, sizeof(char));
 strcpy(temp_ex->name1, name1);
 strcpy(temp_ex->name2, name2);
 strcpy(temp_ex->portname,portname);

 return temp_ex;
}

EXTERNAL *get_external_from_portname(portname, indsys)
char *portname;
SYS *indsys;
{
  EXTERNAL *ext;

  if (strcmp(portname,"") == 0)
    return (EXTERNAL *)NULL;

  for(ext=indsys->externals; ext != (EXTERNAL *)NULL; ext = ext->next)
    if (strcmp(portname,ext->portname) == 0)
      return ext;

  return (EXTERNAL *)NULL;
}

EXTERNAL *get_next_ext(ext)
EXTERNAL *ext;
{
  while(ext != NULL && ext->col_Yindex == -1)
    ext = ext->next;

  return ext;
}

NODES *get_next_treeless_node(node)
NODES *node;
{
  while(node != NULL && (getrealnode(node)->examined != 0))
    node = node->next;

  return node;
}

TREE *make_new_tree()
{
  TREE *temp;

  temp = (TREE *)Gmalloc(sizeof(TREE));
  temp->loops = NULL;
  temp->number_of_loops = 0;
  temp->next = NULL;

  return temp;
}

TREE *add_tree_to_list(tree, t_list)
TREE *tree, *t_list;
{
  tree->next = t_list;
  return tree;
}

NODES *pop_node(stack)
NPATH **stack;
{
  NODES *node;
  NPATH *killme;

  if (stack == NULL)
    return NULL;
  else {
    node = (*stack)->node;
    killme = *stack;
    *stack = (*stack)->next;
    free(killme);

    return node;
  }
}

push_node(node, stack)
NODES *node;
NPATH **stack;
{
  NPATH *newelem;

  newelem = (NPATH *)Gmalloc(sizeof(NPATH));
  newelem->node = node;
  newelem->next = *stack;
  *stack = newelem;
}

/* This function is based on the algorithm from "Graph Theory with Applications
   to Engin. and Comp Sci" by Narsingh Deo. 1974, pp. 280-284. */
make_trees(indsys)
SYS *indsys;
{

  NODES *tempnode, *node, *other, *orig;
  TREE *atree;
  NPATH *stack = NULL;
  SEGLIST *branches;

  indsys->trees = NULL;
  indsys->num_trees = 0;  /* used in fillA() */

  tempnode = get_next_treeless_node(indsys->nodes);
  while(tempnode != NULL) {
    atree = make_new_tree();
    indsys->num_trees++;
    indsys->trees = add_tree_to_list(atree, indsys->trees);
    node = getrealnode(tempnode);
    node->treeptr = atree;
    node->level = 0;
    push_node(node, &stack);
    while(stack != NULL) {
      node = pop_node(&stack);
      /* skip deleted and gp segs, and get original node */
      branches = get_next_branch(node->connected_segs);
      while(branches != NULL) {
	other = getrealnode(getothernode(branches->original, branches->seg));

#if 1==0
   unimplemented junk
	if (!is_extern_seg(branches->seg)) {
	  if (other->level == -1) {
	    other->treeptr = atree;
	    other->level = node->level + 1;
	    other->pred = branches->seg;  /* add branch to tree */
	    push_node(other, &stack);
	  }
	  else {
	    make_loop(node, other, branches->seg, atree); /* make a loop (check for selfloop) */
	  }
	  mark_seg(branches->seg);    /* mark seg as used */
#endif

	if (other->level == -1) {
	  other->treeptr = atree;
	  other->level = node->level + 1;
	  other->pred = branches->seg;  /* add branch to tree */
	  push_node(other, &stack);
	}
	else {
#if 1==0
   junk
	  /* We must treat Voltage sources from .extern statements carefully.
	     This is for the new Preconditioner which requires that a voltage
	     source appear in only one mesh.  The following should keep
	     the pseudo-seg for a source off the tree, but use it only if it
	     completes a mesh.
	  */
	  if (other->level != -1) {
	    make_loop(node, other, branches->seg, atree);
	    mark_seg(branches->seg);
	  }
#endif
	  make_loop(node, other, branches->seg, atree); /* make a loop (check for selfloop) */
	}
#if 1==0
  junk
	branches = get_next_branch(branches->next);

      }  /* end while(branches.. */

#endif
	mark_seg(branches->seg);    /* mark seg as used */
        branches = get_next_branch(branches->next);
      }
      mark_node(node);  /* mark node as examined */
    }
    tempnode = get_next_treeless_node(tempnode);
  }
}

SEGLIST *get_next_branch(b_list)
SEGLIST *b_list;
{
  while(b_list != NULL && (is_marked(b_list->seg) || is_gp(b_list->seg)) )
    b_list = b_list->next;

  return b_list;
}

int is_marked(seg)
seg_ptr seg;
{
  if (seg.type == NORMAL)
    return  ((SEGMENT *)seg.segp)->is_deleted ;
  else if (seg.type == PSEUDO)
    return  ((PSEUDO_SEG *)seg.segp)->is_deleted ;
  else {
    fprintf(stderr, "is_marked: unknown seg type: %d\n",seg.type);
    exit(1);
  }
}

mark_seg(seg)
seg_ptr seg;
{
  if (seg.type == NORMAL)
    ((SEGMENT *)seg.segp)->is_deleted = 1;
  else
    ((PSEUDO_SEG *)seg.segp)->is_deleted = 1;
}

unmark_seg(seg)
seg_ptr seg;
{
  if (seg.type == NORMAL)
    ((SEGMENT *)seg.segp)->is_deleted = 0;
  else
    ((PSEUDO_SEG *)seg.segp)->is_deleted = 0;
}

mark_node(node)
NODES *node;
{
  node->examined = 1;
}

unmark_node(node)
NODES *node;
{
  node->examined = 0;
}

is_node_marked(node)
NODES *node;
{
  return (node->examined != 0);
}

PATHLIST *add_path_to_list(path, list)
SPATH *path;
PATHLIST *list;
{
  PATHLIST *templist;

  templist = (PATHLIST *)Gmalloc(sizeof(PATHLIST));
  templist->path = path;
  templist->next = list;

  return templist;
}

make_loop(node_l, node_s, seg, tree)
NODES *node_s;  /* one branch to main trunk. */
                /* possibly zero branches to main trunk if seg is EXTERNTYPE?
                     (obsolete comment?)*/
NODES *node_l;  /* along main trunk with many branches to where node_s is */
seg_ptr seg;    /* segment connecting above nodes (not in tree) */
TREE *tree;     /* tree that this loop will be contained within */
{

  SPATH *path = NULL;
  seg_ptr pre_seg;
  NODES *pre_node;
  int count;

  path = add_seg_to_list(seg, path); /* make a path */

  if (node_s != node_l) { /* is this not a self loop? */
    pre_node = node_l;
    count = node_l->level - (node_s->level - 1);
#if 1==0
   junk
    if (is_extern_seg(seg))
      count--;  /* let's assume node_s is on main trunk, */
                /* so this makes up for -1 above*/
#endif
    if (count < 0) {
      fprintf(stderr, "make_loop:  strange loop. count < 0?\n");
      exit(1);
    }
    while(count > 0) {
      pre_seg = pre_node->pred;
      path = add_seg_to_list(pre_seg, path);
      pre_node = getrealnode(getothernode(pre_node, pre_seg));
      count--;
    }
#if 1==0
   junk
    if (!is_extern_seg(seg)) {
      pre_seg = node_s->pred;
      path = add_seg_to_list(pre_seg, path);
      if (getrealnode(getothernode(node_s, pre_seg)) != pre_node) {
	fprintf(stderr, "Internal Error in make_loop: Hey, these don't make a loop!\n");
	exit(1);
      }
#endif
    pre_seg = node_s->pred;
    path = add_seg_to_list(pre_seg, path);
    if (getrealnode(getothernode(node_s, pre_seg)) != pre_node) {
      fprintf(stderr, "make_loop: Hey, these don't make a loop!\n");
      exit(1);
    }
#if 1==0
  junk
    else if (pre_node != node_s) {
      if (pre_node->level == node_s->level &&
	  getrealnode(getothernode(node_s, node_s->pred))
	  == getrealnode(getothernode(pre_node, pre_node->pred)) ) {
	/* the assumption that node_s was on the main trunk was wrong */
	/*  We need two more segments (note: order is important)*/
	path = add_seg_to_list(pre_node->pred, path);
	path = add_seg_to_list(node_s->pred, path);
      }
      else {
	fprintf(stderr, "Internal Error in make_loop: Hey, these don't make a loop (extern)!\n");
	exit(1);
      }
    }
#endif
  }

  tree->loops = add_path_to_list(path, tree->loops);
  tree->number_of_loops += 1;
  /* add this to other things too? */
}

int count_tree_meshes(trees)
TREE *trees;
{
  int total = 0;

  while(trees != NULL) {
    total += trees->number_of_loops;
    trees = trees->next;
  }

  return total;
}

#if 1==0
  unimplemented obsolete junk

/* This estimates the number of extra meshes to be produced by breaking
   all the big meshes into many smaller ones which will have at most
   fils_per_mesh filaments per mesh
*/
estimate_extra_meshes(trees, fils_per_mesh)
TREE *trees;
int fils_per_mesh;
{
  int total = 0;
  PATHLIST *plist;
  SPATH *one_seg;
  int segs_in_loop;

  while(trees != NULL) {
    for(plist = trees->loops; plist != NULL; plist = plist->next) {
      segs_in_loop = 0;
      for(one_seg = plist->path; one_seg != NULL; one_seg = one_seg->next) {
	if (is_normal_seg(one_seg->seg))
	  segs_in_loop += fils_per_mesh; /*a gross overestimate. 1 is better*/
	else if (is_pseudo_seg(one_seg->seg))
	  segs_in_loop += get_upper_num_segs((PSEUDO_SEG *)one_seg->seg.segp);
	else {
	  fprintf(stderr,"Internal Err: bad seg type in estimate_ext...\n");
	  exit(1);
	}
      }
      total += segs_in_loop/fils_per_mesh + 1;
    }
    trees = trees->next;
  }
  return total;
}
#endif

int count_externals(ext_list)
EXTERNAL *ext_list;
{
  int count = 0;

  while(ext_list != NULL) {
    ext_list = ext_list->next;
    count++;
  }

  return count;
}

#if 1==0
get_upper_num_segs(pseg)
PSEUDO_SEG *pseg;
{
  return pseg->upper_num_segs;
}
#endif //1==0, Enrico


/****  The following are functions for finding meshes created by holes ***/
/****  in ground planes.  Many functions are near duplicates of those above **/

/* this finds all the extra meshes that result from holes in the plane.
   Each hole needs one mesh which outlines it.
   It adds a tree to indsys->trees which contains all the new loops.
*/

find_hole_meshes(indsys)
SYS *indsys;
{
  NODES *node, *tnode, ***pnodes;
  TREE *atree;
  NPATH *nodes_in_hole, *surrounding_nodes;
  NPATH *stack = NULL;
  GROUNDPLANE *plane;
  int i, j, nodes1, nodes2, s1, s2;

  atree = make_new_tree();

  node = get_next_gphole_node(indsys->nodes);
  nodes_in_hole = NULL;

  while(node != NULL) {
    push_node(node, &stack);
    nodes_in_hole = add_node_to_list(node, nodes_in_hole);
    mark_node(node);
    while(stack != NULL) {
      node = pop_node(&stack);
      s1 = node->s1;
      s2 = node->s2;
      plane = node->gp;
      pnodes = plane->pnodes;
      nodes1 = plane->num_nodes1;
      nodes2 = plane->num_nodes2;

      /* check node's eight neighbors for hole nodes and add to list */
      for(i = MAX(s1 - 1,0); i <= MIN(s1 + 1, nodes1 - 1); i++)
	for(j = MAX(s2 - 1,0); j <= MIN(s2 + 1, nodes2 - 1); j++)
	  if (!(i == s1 && j == s2)) {
	    tnode = pnodes[i][j];
	    if (is_hole(tnode) && !is_node_marked(tnode)) {
	      push_node(tnode, &stack);
	      nodes_in_hole = add_node_to_list(tnode, nodes_in_hole);
	      mark_node(tnode);
	    }
	  }
    }

    /* find all surrounding nodes */
    surrounding_nodes = find_surrounding(nodes_in_hole);

    clear_marks_and_level(surrounding_nodes);

    /* find all the circuits formed.  Clear seg->is_deleted also */
    make_gp_trees(surrounding_nodes, atree);

    free_nodelist(nodes_in_hole);
    free_nodelist(surrounding_nodes);
    nodes_in_hole = NULL;
    node = get_next_gphole_node(node);
  }

  if (atree->number_of_loops != 0) {
    indsys->trees = add_tree_to_list(atree, indsys->trees);
    mark_used_segs(atree->loops);
  }
  else
    free(atree);

}

NODES *get_next_gphole_node(node)
NODES *node;
{
  while(node != NULL && (!is_hole(node) || is_node_marked(node)))
    node = node->next;

  return node;
}

NPATH *find_surrounding(nodes_in_hole)
NPATH *nodes_in_hole;
{
  NPATH *np;
  NPATH *surround = NULL;
  int s1, s2, nodes1, nodes2, i, j;
  NODES ***pnodes;
  NODES *node, *tnode;
  GROUNDPLANE *plane;

  for(np = nodes_in_hole; np != NULL; np = np->next) {
    node = np->node;
    s1 = node->s1;
    s2 = node->s2;
    plane = node->gp;
    pnodes = plane->pnodes;
    nodes1 = plane->num_nodes1;
    nodes2 = plane->num_nodes2;

    /* check node's eight neighbors for hole nodes */
    for(i = MAX(s1 - 1,0); i <= MIN(s1 + 1, nodes1 - 1); i++)
      for(j = MAX(s2 - 1,0); j <= MIN(s2 + 1, nodes2 - 1); j++)
	if (!(i == s1 && j == s2)) {
	  tnode = pnodes[i][j];
	  if (!is_hole(tnode) && !is_orignode_in_list(tnode, surround))
	    surround = add_node_to_list(tnode, surround);
	}
  }

  return surround;
}

clear_marks_and_level(nlist)
NPATH *nlist;
{
  NPATH *np;

  for(np = nlist; np != NULL; np = np->next) {
    unmark_node(np->node);
    np->node->level = -1;
  }
}

/* find all the circuits (loops) formed by the nodes in nlist.
   This is nearly identical to make_trees().
   Also clear seg->is_deleted when done */
make_gp_trees(nlist, atree)
NPATH *nlist;
TREE *atree;
{
  NPATH *np;
  NODES *node, *other;
  NPATH *stack = NULL;
  SPATH *path;
  int connected_counter;
  seg_ptr seg;
  int counter = 0;

  seg.type = NORMAL;

  np = get_next_unexamined_node(nlist);
  while(np != NULL) {
    node = np->node;
    node->level = 0;
    push_node(node, &stack);
    while(stack != NULL) {
      node = pop_node(&stack);
      connected_counter = 0;
      seg.segp = (void *)get_next_around_hole(node, &connected_counter, nlist);
      while(seg.segp != NULL) {
	other = getothernode(node, seg);
	if (other->level == -1) {
	  other->level = node->level + 1;
	  other->pred = seg;
	  push_node(other, &stack);
	}
	else {
	  path = make_gp_loop(node, other, seg);
	  counter++;
	  atree->loops = add_path_to_list(path, atree->loops);
	  atree->number_of_loops += 1;
	}
	mark_seg(seg);
	seg.segp = (void *)get_next_around_hole(node,&connected_counter,nlist);
      }
      mark_node(node);
    }
    np = get_next_unexamined_node(np);
  }

  if (counter > 1) {
    printf("Warning: Multiple boundaries found around one hole region\n");
    printf("  possibly due to an isolated or nearly isolated region of conductor.\n");
    printf("  This may lead to no unique solution.\n");
  }

  /* we'll clear all the loops due to ground planes each time since
     hopefully there won't be many */
  clear_used_segs(atree->loops);
}

NPATH *get_next_unexamined_node(np)
NPATH *np;
{
  while(np != NULL && is_node_marked(np->node))
    np = np->next;

  return np;
}

/* this gets the next segment which is not marked and is surrounding the hole*/
SEGMENT *get_next_around_hole(node, counter, nlist)
NODES *node;
int *counter;
NPATH *nlist;
{
  static seg_ptr seg;

  seg.segp = (void *)get_next_gp_seg(node, counter);

  while(seg.segp != NULL
	&& (is_marked(seg)
	    || !is_orignode_in_list(getothernode(node,seg),nlist)))
      seg.segp = (void *)get_next_gp_seg(node, counter);

  return (SEGMENT *)seg.segp;
}

/* returns each of the four segments around a ground plane node (if they
   exist). */
SEGMENT *get_next_gp_seg(node, counter)
NODES *node;
int *counter;
{
  GROUNDPLANE *plane = node->gp;
  int s1 = node->s1;
  int s2 = node->s2;

  (*counter)++;

  /* start with seg to the right */
  if (*counter == 1) {
    if (s1 < plane->num_nodes1 - 1 && plane->segs1[s1][s2] != NULL)
      return plane->segs1[s1][s2];
    else
      (*counter)++;
  }

  /* next is seg above */
  if (*counter == 2) {
    if (s2 < plane->num_nodes2 - 1 && plane->segs2[s1][s2] != NULL)
      return plane->segs2[s1][s2];
    else
      (*counter)++;
  }

  /* next is to the left */
  if (*counter == 3) {
    if (s1 > 0 && plane->segs1[s1-1][s2] != NULL)
      return plane->segs1[s1-1][s2];
    else
      (*counter)++;
  }

  /* next is below */
  if (*counter == 4) {
    if (s2 > 0 && plane->segs2[s1][s2-1] != NULL)
      return plane->segs2[s1][s2-1];
    else
      (*counter)++;
  }

  return NULL;

}

/* given a starting and ending node, this makes a path from the tree
   information.  This is almost identical to make_loops().  I just had
   to take out all the 'getrealnode' calls and I made it return a SPATH. */
SPATH *make_gp_loop(node_l, node_s, seg)
NODES *node_s;  /* one branch to main trunk */
NODES *node_l;  /* along main trunk with many branches to where node_s is */
seg_ptr seg;    /* segment connecting above nodes (not in tree) */
{

  SPATH *path = NULL;
  seg_ptr pre_seg;
  NODES *pre_node;
  int count;

  path = add_seg_to_list(seg, path); /* make a path */

  if (node_s != node_l) { /* is this not a self loop? */
    pre_node = node_l;
    count = node_l->level - (node_s->level - 1);
    if (count < 0) {
      fprintf(stderr, "make_loop:  strange loop. count < 0?\n");
      exit(1);
    }
    while(count > 0) {
      pre_seg = pre_node->pred;
      path = add_seg_to_list(pre_seg, path);
      pre_node = getothernode(pre_node, pre_seg);
      count--;
    }
    pre_seg = node_s->pred;
    path = add_seg_to_list(pre_seg, path);
    if (getothernode(node_s, pre_seg) != pre_node) {
      fprintf(stderr, "make_loop: Hey, these don't make a loop!\n");
      exit(1);
    }
  }

  return path;
}

clear_used_segs(plist)
PATHLIST *plist;
{
  PATHLIST *pl;
  SPATH *path;

  for(pl = plist; pl != NULL; pl = pl->next)
    for(path = pl->path; path != NULL; path = path->next)
      unmark_seg(path->seg);
}

/*this marks used segs so that path_through_gp() has a record of these meshes*/
mark_used_segs(plist)
PATHLIST *plist;
{
  PATHLIST *pl;
  SPATH *path;

  for(pl = plist; pl != NULL; pl = pl->next)
    for(path = pl->path; path != NULL; path = path->next)
      ((SEGMENT *)path->seg.segp)->is_deleted++;
}


/* find the nearest node to node in nodelist */
NODES *find_nearest_node(node, nodelist)
NODES *node;
NPATH *nodelist;
{
  NODES *min_node;
  double min, dist;

  if (nodelist == NULL)
    return NULL;

  min_node = nodelist->node;
  min = dist_between_nodes(node, min_node);

  nodelist = nodelist->next;

  while(nodelist != NULL) {
    dist = dist_between_nodes(node, nodelist->node);
    if (dist < min) {
      min = dist;
      min_node = nodelist->node;
    }
    nodelist = nodelist->next;
  }

  return min_node;
}

double dist_between_nodes(n1, n2)
NODES *n1, *n2;
{
  return magsq( n1->x - n2->x, n1->y - n2->y, n1->z - n2->z);
}





