/*
 * $Id$
 */

/*
 * Binary Tree
 *  ~~~~~~ ~~~~
 * Very general purpose routines for binary tree implemetation.
 * Each leaf contains a (char*)key with the name of the leaf
 * and a (void*)value which contains the value of the leaf.
 *
 * The searching is done with the key's checked with strcmp
 * that means that an INTEGER or a REAL is stored according
 * to its binary representation in memory.
 *
 * When adding a leaf no memory allocation is done for the key
 * and the value.
 */

#include <ostream>
#include <string.h>

#include "os.h"
#include "bintree.h"

using namespace std;

#ifdef _DEBUG
//static int ScanDepth(BinLeaf *leaf, int depth);
#endif

/** BinTree
 * Copy one binary tree to another, do not create new
 * data, just points the other tree!!!
 *
 * WARNING:
 * 1. This can generate problems if someone decides
 *    to delete the data!!!
 * 2. The destination tree must be empty!!!
 */
BinTree::BinTree(BinTree& src)
{
	copyLeaf(src.parent);
} // BinTree

/** add new entry
 * @param name	pointer to key. IMPORTANT the bintree doesn't own the key
 * @param dat	item to be stored
 * @return BinLeaf position
 */
BinLeaf* BinTree::add(const char *name, void *dat)
{
	BinLeaf *thisEntry;
	BinLeaf *lastEntry=NULL;
	BinLeaf *leaf;
	int	leftTaken=0;
	int	dep=0;

	/* If tree is NULL then it will produce an error */
	thisEntry = parent;
	while (thisEntry != NULL) {
		lastEntry = thisEntry;
		int cmp = STRCMP(name,thisEntry->key);
		if (cmp < 0) {
			thisEntry = thisEntry->left;
			leftTaken = true;
		} else
		if (cmp > 0) {
			thisEntry = thisEntry->right;
			leftTaken = false;
		} else
			return thisEntry;
		dep++;
	}

	/* Create a new entry */
	leaf = new BinLeaf(name,dat);

	if (parent==NULL)
		parent = leaf;
	else {
		if (leftTaken)
			lastEntry->left = leaf;
		else
			lastEntry->right = leaf;
	}
	nitems++;
	if (dep>maxDepth) {
		maxDepth = dep;
		if (maxDepth > balanceDepth)
			balance();
	}
	return leaf;
} // add

/** find
 * @param name	key to search
 * @return value of key or NULL
 */
void* BinTree::find(const char *name) const
{
	BinLeaf *leaf;

	leaf = parent;
	while (leaf != NULL) {
		int cmp = STRCMP(name, leaf->key);
		if (cmp < 0)
			leaf = leaf->left;
		else
		if (cmp > 0)
			leaf = leaf->right;
		else
			return leaf->data;
	}
	return NULL;
} // Find

/** delete an item from the tree
 * @param name	key to be removed.
 * @return true/false if the key was found
 * @warning if DelFunc is set it will call
 * it to free the memory otherwise it does nothing
 */
/* -------------------------------------------------------------- */
/* To correctly delete a pointer from a Binary tree we must not   */
/* change the tree structure, that is the smaller values are the  */
/* left most. In order to satisfy this with few steps we must     */
/* replace the pointer that is to be erased with the one which is */
/* closest with a smaller value (the right most from the left     */
/* branch, as you can see below                                   */
/*                     ...                                        */
/*                    /                                           */
/*                 (name)   <-- to be dropped                     */
/*                 /    \                                         */
/*              (a)      (d)       (c)= newid from left branch    */
/*             /   \       ...          where  c->right=NULL      */
/*          ...     (c)            (a)= par_newidt parent of newid*/
/*                 /   \                                          */
/*              (b)     NIL       newid will become the new       */
/*             ...                sub-head when (name) is dropped */
/*                    |                                           */
/*                    |                                           */
/*                   \|/                                          */
/*                    V   ...                                     */
/*                      /                                         */
/*                    (c)           but in the case that          */
/*                   /   \          (a)=(c) when a->right = NULL  */
/*                 (a)    (d)       then the tree is very simple  */
/*                /   \    ....     we simply replace a->right=d  */
/*              ...   (b)                                         */
/*                    ....                                        */
/*                                                                */
/* -------------------------------------------------------------- */
bool BinTree::del(const char *name)
{
	BinLeaf *thisid,
		 *previous = NULL,
		 *par_newid,
		 *newid;
	bool	leftTaken=0;

	thisid = parent;
	while (thisid != NULL) {
		int cmp = STRCMP(name,thisid->key);
		if (cmp < 0) {
			previous  = thisid;
			thisid    = thisid->left;
			leftTaken = true;
		} else
		if (cmp > 0) {
			previous  = thisid;
			thisid    = thisid->right;
			leftTaken = false;
		} else
			break;
	}
	if (thisid == NULL) return false;	/* Not Found */

	if (thisid->right == NULL)
		newid = thisid->left;
	else
	if (thisid->left == NULL)
		newid = thisid->right;
	else {		/* when no node is empty */
		/* find the right most id of the */
		/* left branch of thisid         */

		par_newid = thisid;
		newid     = thisid->left;
		while (newid->right != NULL) {
			par_newid = newid;
			newid = newid->right;
		}

		/* newid must now replace thisid */
		newid->right = thisid->right;

		if (par_newid != thisid) {
			par_newid->right = newid->left;
			newid->left = thisid->left;
		}
	}

	if (thisid == parent)
		parent = newid;
	else {
		if (leftTaken)
			previous->left = newid;
		else
			previous->right = newid;
	}
	thisid->left = NULL;
	thisid->right = NULL;
	destroyBranch(thisid);

	return true;
} // del

/** destroyBranch  */
void BinTree::destroyBranch(BinLeaf *leaf)
{
	if (!leaf) return;

	if (leaf->left)
		destroyBranch(leaf->left);

	if (leaf->right)
		destroyBranch(leaf->right);

	if (leaf->data && _delFunc!=NULL)
		_delFunc(leaf->key, leaf->data);

	if (leaf==parent)
		parent = NULL;

	delete leaf;
	nitems--;
} // destroyBranch

/** copyLeaf */
void BinTree::copyLeaf(BinLeaf* src)
{
	if (src == NULL)
		return;

	/* add them in a sorted way */
	copyLeaf(src->left);
	add(src->key, src->data);
	copyLeaf(src->right);
} // copyLeaf

/* forEachLeaf
 * stop if the traverse function returns true
 */
int BinTree::forEachLeaf(BinLeaf *leaf, BinTraverse func, void *data)
{
	if (leaf == NULL)
		return false;

	if (forEachLeaf(leaf->left,func,data))
		return true;

	if ((func)(leaf->key,leaf->data,data))
		return true;

	return forEachLeaf(leaf->right,func,data);
} // forEachLeaf

/** scanDepthLeaf */
int BinTree::scanDepthLeaf(BinLeaf *leaf, int depth)
{
	int left, right;
	if (leaf==NULL)
		return depth;
	left  = scanDepthLeaf(leaf->left,depth+1);
	right =	scanDepthLeaf(leaf->right,depth+1);
	return MAX(left,right);
} // scanDepthLeaf

/** printLeaf */
void BinTree::printLeaf(ostream& os, BinLeaf *leaf, int depth) const
{
	long i;

	if (!leaf) return;

	printLeaf(os,leaf->left,depth+3);

	for (i=0; i<depth; i++) os<< ' ';
	os << '\"' << leaf->key << '\"';

//	if (leaf->data)
//		printf("%p\n",leaf->data);
//	else
//		printf("NULL\n");

	printLeaf(os,leaf->right,depth+3);
} // printLeaf

/** balanceLeaf */
void BinTree::balanceLeaf(BinLeaf *leaf, BinLeaf **head, BinLeaf **tail)
{
	BinLeaf *Lhead, *Ltail;
	BinLeaf *Rhead, *Rtail;

	if (leaf == NULL) {
		*head = NULL;
		*tail = NULL;
		return;
	}

	balanceLeaf(leaf->left,  &Lhead, &Ltail);
	balanceLeaf(leaf->right, &Rhead, &Rtail);

	/* connect nodes */
	/*  head - left - middle - right - tail */

	if (Ltail) Ltail->right = leaf;
	leaf->left   = Ltail;

	leaf->right  = Rhead;
	if (Rhead) Rhead->left  = leaf;

	if (Lhead)
		*head = Lhead;
	else
		*head = leaf;

	if (Rtail)
		*tail = Rtail;
	else
		*tail = leaf;
} // balanceLeaf

/** constructLeaf */
BinLeaf *BinTree::constructLeaf(BinLeaf *head, BinLeaf *tail, int n, int *maxd)
{
	int	Lmaxd, Rmaxd, i, mid;
	BinLeaf *Lleaf,
		 *Rleaf,
		 *LMidleaf,
		 *Midleaf,
		 *RMidleaf;

	if (n==0) return NULL;
	if (n==1) {
		/* then head must be equal to tail */
		head->left  = NULL;
		head->right = NULL;
		return head;
	}
	if (n==2) {
		(*maxd)++;
		head->left  = NULL;
		head->right = tail;
		tail->left  = NULL;
		tail->right = NULL;
		return head;
	}

	/* --- find middle --- */
	mid = n/2;
	LMidleaf = head;
	for (i=0; i<mid-1; i++)
		LMidleaf = LMidleaf->right;
	Midleaf = LMidleaf->right;
	RMidleaf = Midleaf->right;

	/* --- do the same for left and right branch --- */
	Lmaxd = Rmaxd = *maxd+1;

	Lleaf = constructLeaf(head,LMidleaf,mid,&Lmaxd);
	Rleaf = constructLeaf(RMidleaf,tail,n-mid-1,&Rmaxd);

	*maxd = MAX(Lmaxd, Rmaxd);

	Midleaf->left = Lleaf;
	Midleaf->right = Rleaf;

	return Midleaf;
} // constructLeaf

/** balance the binary tree */
void BinTree::balance()
{
	BinLeaf *head,
		 *tail;
	int	  maxd=1;

	/* make tree a double queue */
	balanceLeaf( parent, &head, &tail );

	/* reconstruct the tree */
	parent       = constructLeaf( head, tail, nitems, &maxd );
	maxDepth     = maxd;
	balanceDepth = maxDepth+BALANCE_INC;
} // balance
