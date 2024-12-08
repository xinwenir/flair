/*
 * $Id$
 */

#ifndef __BINTREE_H
#define __BINTREE_H

#define BALANCE_START	7
#define BALANCE_INC	8

#include <stdlib.h>
#include <ostream>

typedef int  (BinTraverse)(const char *key, void *value, void *data);
typedef void (BinDel)(const char *key, void *value);

/** Binary Tree Leaf data structure */
class BinLeaf {
protected:
const	char	*key;
	void	*data;
	BinLeaf	*left,
		*right;
public:
	BinLeaf(const char *n, void *d) : key(n), data(d) {
		left = NULL;
		right = NULL;
	}
const	char*	getKey()	const { return key; }
	void*	getData()	const { return data; }
	void*	getLeft()	const { return left; }
	void*	getRight()	const { return right; }
	friend class BinTree;
}; // class BinLeaf

/** Binary Tree */
class BinTree {
private:
	BinLeaf	*parent;
	int	 nitems;
	int	 maxDepth;
	int	 balanceDepth;
	BinDel	*_delFunc;

public:
	BinTree() {
		parent       = NULL;
		nitems       = 0;
		maxDepth     = 0;
		_delFunc     = NULL;
		balanceDepth = BALANCE_START;
	}
	BinTree(BinTree& src);

	BinTree(BinDel func) : _delFunc(func) {
		parent       = NULL;
		nitems       = 0;
		maxDepth     = 0;
		balanceDepth = BALANCE_START;
	}

	~BinTree()		{ destroyTree(); }

	BinLeaf* add(const char *name, void *dat);
	void*	find(const char *name) const;
	/** @return item with key name */
	void*	operator [](const char *name) const { return find(name); }
	bool	has(const char *name) { return find(name)!=NULL; }
	bool	del(const char *name);

	/** @return number of items in tree */
	int	count() const { return nitems; }

	void	balance();

	/** destroy the whole tree */
	void	destroyTree() { destroyBranch(parent); }

	/** traverse every item of tree
	 * @param func	function to call
	 * @param data	extra data to be passed to func
	 */
	void	forEach(BinTraverse func,void *data)
			{ forEachLeaf(parent,func,data); }

	/** @return the maximum depth of the tree */
	int	scanDepth()
			{ return scanDepthLeaf(parent,0); }

	/** print tree on a stream */
	void	printOn(std::ostream& os) const
			{ printLeaf(os,parent,0); }

	/** set function to call when a key is deleted */
	void	deleteFunc(BinDel func)		{ _delFunc = func; }

	/** @return memory used by bintree (only of the structrure) */
	size_t	memory() const	{ return sizeof(BinTree) + count()*sizeof(BinLeaf); }

private:
	void	copyLeaf(BinLeaf* src);
	int	forEachLeaf(BinLeaf* leaf, BinTraverse func, void *data);
	int	scanDepthLeaf(BinLeaf *leaf, int depth);
	void	balanceLeaf(BinLeaf *leaf, BinLeaf **head, BinLeaf **tail);
	BinLeaf *constructLeaf(BinLeaf *head, BinLeaf *tail, int n, int *maxDepth);
	void	printLeaf(std::ostream& os, BinLeaf *leaf, int depth) const;
	void	destroyBranch(BinLeaf *leaf);
}; // class BinTree

// Output to a stream.
inline std::ostream& operator << (std::ostream& os, const BinTree& tree)
{
	tree.printOn(os);
	return os;
}
#endif
