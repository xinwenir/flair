/*
 * $Id$
 *
 * Copyright and User License
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright 2006-2019 CERN and INFN
 * 
 *
 * Please consult the LICENSE file for the license 
 *
 * DISCLAIMER
 * ~~~~~~~~~~
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
 * NOT LIMITED TO, IMPLIED WARRANTIES OF MERCHANTABILITY, OF
 * SATISFACTORY QUALITY, AND FITNESS FOR A PARTICULAR PURPOSE
 * OR USE ARE DISCLAIMED. THE COPYRIGHT HOLDERS AND THE
 * AUTHORS MAKE NO REPRESENTATION THAT THE SOFTWARE AND
 * MODIFICATIONS THEREOF, WILL NOT INFRINGE ANY PATENT,
 * COPYRIGHT, TRADE SECRET OR OTHER PROPRIETARY RIGHT.
 *
 * LIMITATION OF LIABILITY
 * ~~~~~~~~~~~~~~~~~~~~~~~
 * THE COPYRIGHT HOLDERS AND THE AUTHORS SHALL HAVE NO
 * LIABILITY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
 * CONSEQUENTIAL, EXEMPLARY, OR PUNITIVE DAMAGES OF ANY
 * CHARACTER INCLUDING, WITHOUT LIMITATION, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES, LOSS OF USE, DATA OR PROFITS,
 * OR BUSINESS INTERRUPTION, HOWEVER CAUSED AND ON ANY THEORY
 * OF CONTRACT, WARRANTY, TORT (INCLUDING NEGLIGENCE), PRODUCT
 * LIABILITY OR OTHERWISE, ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGES.
 *
 * Author:	Vasilis.Vlachoudis@cern.ch
 * Date:	04-Feb-2002
 */

#ifndef __ARRAY_H
#define __ARRAY_H

#include <assert.h>
#include <stdlib.h>
#include <string.h>

template <typename T>
struct _is_pointer {
	static const bool value = false;
};

template <typename T>
struct _is_pointer<T*> {
	static const bool value = true;
};

// C++11 makes it obsolete it has a std::is_pointer()
template <typename T>
bool IsPointer(const T&) {
	return _is_pointer<T>::value;
}

class ArrayEmpty {};

/*-------------------------------------------- */
/** A dynamic array sorted or unsorted of pointers to objects.
 * The array does not own the contents of the data.
 * @see ArrayIterator
 */
template <class T>
class Array
{
	/** Compare function for sorted arrays must return -1,0,+1	*/
	typedef int (*ArrayCompare)(const T&, const T&);

protected:
	T		*data;		/** data array			*/
	int		_capacity;	/** size of allocated array	*/
	int		_count;		/** number of data in array	*/
	int		_delta;		/** increment of size		*/
	ArrayCompare	_compare;	/** compare function		*/

public:
	/** Initialize the array with a predefined size
	 * @param sz initial size
	 * @param d=0 delta step
	 */
	Array(int d=8)
		: data(NULL),
		  _capacity(0),
		  _count(0),
		  _delta(d),
		  _compare(NULL)
		{ }

	/** Copy constructor */
	Array(const Array<T>&);

	/** Delete array but not data */
	~Array() { if (data) delete [] data; }

	/** Copy array */
	Array<T>& operator = (const Array<T>&);

	/**
	 * @param index index of data to be returned
	 * @return data pointer */
	T&	operator [] (int index) const {
			assert(_capacity > 0 && data != NULL && index < _count);
			return data[index];
		}

	/**
	 * @param index index of data to be returned
	 * @return data pointer */
	T&	get(int index) const {
			assert(_capacity > 0 && data != NULL && index < _count);
			return data[index];
		}

	/** allocate a fixed array size
	 * @param n	size of array
	 */
	void	allocate(int n)	{
			resize(n);
			forceCount(n);
		}

	/** force count to argument. WARNING very dangerous */
	void	forceCount(int c) {
			assert(0<=c && c<=_capacity);
			_count = c;
		}

	/** resize the array starting at offset
	 * @param newSz the new size of the array
	 */
	bool	resize(int newSz);

	/** Define a compare function to create a sorted list.
	 * The function should be set before populating any data.
	 * Otherwise use the sort method to sort the data.
	 * @param cmp function that compares two data pointers
	 *		should return -1, 0 or 1
	 * @see sort to sort already stored data
	 */
	void	compare(ArrayCompare cmp=NULL) { _compare = cmp; }

	ArrayCompare compare() { return _compare; }

	/** copy the data from another pointer */
	void	copyData(T* src, int n);

	/** shift data from position idx with length len
	 * @param idx location
	 * @param len to shift
	 */
	void	shiftby(int idx, int len);

	/** add a new data pointer t at the end of the array. If
	 * the array is sorted add it to the correct position.
	 * @param t data pointer to add
	 * @return position
	 */
	int	add(const T& t);
	int	append(const T& t)	{ return add(t); }
	int	push_back(const T& t)	{ return add(t); }

	Array<T>& operator <<(const T& t)	{ add(t); return *this; }

	/** insert a new data pointer at specific location
	 * @param idx location
	 * @param t data pointer to add
	 */
	void	insert(int idx, const T &t) {
			shiftby(idx,1);
			data[idx] = t;
		}

	/** equivalent to add
	 * @param t pointer to add
	 * @return position of addition
	 * @see add */
	int	push(const T &t)	{ return add(t); }

	/** removes the last argument */
//	 * @return the last argument */
//	T&	pop(void) {
	void	pop(void) {
			if (!_count) throw ArrayEmpty();
			detach(top());
		}

	/** @return first element */
const	T&	head(void) const {
			if (_count>0) return data[0];
			else throw ArrayEmpty();
		}

	/** @return first element */
	T&	head(void) {
			if (_count>0) return data[0];
			else throw ArrayEmpty();
		}

	/** @return last element const */
const	T&	tail(void) const {
			if (_count>0) return data[_count-1];
			else throw ArrayEmpty();
		}

	/** @return last element */
	T&	tail(void) {
			if (_count>0) return data[_count-1];
			else throw ArrayEmpty();
		}

	/** replace data pointer with a new one
	 * @param idx location where to add data
	 * @param t new data to replace at */
	void	replace(int idx, const T &t);

	/** find and remove data t
	 * @param t data to be removed */
	void	detach(const T &t)	{ erase(t); }
	void	erase(const T &t);

	/** remove data from location
	 * @param idx location of data to be removed */
	void	detach(int idx)		{ erase(idx); }
	void	erase(int idx);
	void	erase(int from, int to);

	/** remove all data */
	void	detachAll()	{ clear(); }
	void	clear();

	/** bubble sort array if a compare function is set. This
	 * method should be called after any manipulation
	 * of the array that will change the order of the data.
	 * @see compare method
	 */
	void	sort();

	/** find position of the data t
	 * if compare function exists it will find data by content
	 * @param t data to search for
	 * @return location of data if succeed
	 *	   -1 if smaller than [0] element
	 *	   count() if bigger than [top()] element
	 * @see compare
	 */
	int	find(const T &t) const;

	/** search for the closest index (lower) when compare is present
	 * @param t data to search for
	 * @return lower bound of data
	 * @see compare
	 */
	int	search(const T &t) const;

	/** find position of the data-pointer t
	 * @param t data to search for
	 * @return location of data
	 */
	int	findPtr(const T &t) const;

	/** @return true if array is full */
	bool	full()	const { return ((_count == _capacity) && _delta == 0); }

	/** @return true if array is empty */
	int	empty()	const { return _count == 0; }

	/** @return the allocated capacity of the array */
	int	capacity() const { return _capacity; }

	/** @return the number of items in the array */
	int	count()	const { return _count; }
	int	size()	const { return _count; }

	/** @return the topmost item # of the list */
	int	top()	const { return _count-1; }

	/** @return the increase delta step of the array */
	int	delta()	const { return _delta; }

	/** set delta increment value
	 * @param d increment step
	 */
	void	delta(const int d) { _delta=d; }

	/** move an existing element to the end of the array */
	void    moveToEnd(size_t idx);
	void    moveToEnd(T element);

	/** @return memory used by structure */
	size_t	memory() const	{ return sizeof(Array<T>) + capacity()*sizeof(T); }

protected:
	/** set to NULL a range [lwr, upr] in the array
	 * @param lwr lower limit
	 * @param upr upper limit */
	void	zero(int lwr, int upr);
	void	zero(int idx);

}; // Array

/*------------------------------------------------------------------------*/
/** Iterator class for the Array class
 * @see Array
 */
template<class T>
class ArrayIterator
{
private:
const	Array<T> *vect;
	int	cur;
	int	lower, upper;

public:
	/** Iterator constructor
	 * @param v array to scan
	 */
	ArrayIterator(const Array<T>& v) {
		vect = &v;
		restart(0, v.count());
	}

	/** Iterator constructor
	 * @param v array to scan
	 * @param start starting element
	 * @param stop end element
	 */
	ArrayIterator(	const Array<T>& v,
			int start,
			int stop)
	{
		vect = &v;
		restart(start, stop);
	}

	/** Initialize iterator
	 * @param v array to scan
	 * @param start starting element
	 * @param stop end element
	 */
	void init(	const Array<T>& v,
			int start=0,
			int stop=0)
	{
		vect = &v;
		if (stop==0) stop = v.count();
		restart(start, stop);
	}

	/** Test if end has reached
	 * @return true if the end has not been reached
	 */
	operator int() { return cur < upper; }

	/** @return current data pointer */
	T &current()
		{ return (cur < upper) ? (*vect)[cur] : (*vect)[upper-1]; }

	/** @return current index */
	int index() const { return cur; }

	/** Increase by one the iterator
	 * @return new data pointer
	 */
	T &operator ++ (int) {
		if (cur >= upper)
			return (*vect)[upper-1];
		else
			return (*vect)[cur++];
	}

	/** Increase by one the iterator
	 * @return new data pointer
	 */
	T &operator ++ () {
		if (cur < upper)
			cur++;
		if (cur >= upper)
			return (*vect)[upper-1];
		else
			return (*vect)[cur];
	}

	/** restart scanning */
	void restart() { restart(lower, upper); }

	/** restart scanning in a specific range
	 * @param start starting element
	 * @param stop ending element
	 */
	void restart(int start, int stop) {
		cur = lower = start;
		upper = stop;
	}
}; // ArrayIterator

/*------------------------------------------------------------------------*/
/** Inverse Iterator class for the Array class
 * @see Array
 */
template<class T>
class InverseArrayIterator
{
private:
const	Array<T>	*vect;
	int	cur;
	int	lower, upper;

public:
	/** Iterator constructor
	 * @param v array to scan
	 */
	InverseArrayIterator(const Array<T>& v) {
		vect  = &v;
		upper = v.count()-1;
		lower = 0;
		cur   = upper;
	}

	InverseArrayIterator(const Array<T>& v, int start, int stop) {
		vect  = &v;
		upper = start;
		lower = stop;
		cur   = upper;
	}

	/** Test if end has reached
	 * @return true if the end has not been reached
	 */
	operator int() { return cur >= lower; }

	/** @return current data pointer */
	T &current() { return (*vect)[cur]; }

	/** @return current index */
	int index() const { return cur; }

	/** Increase by one the iterator
	 * @return new data pointer
	 */
	T &operator ++ (int) {
		assert(cur >= lower);
		return (*vect)[cur--];
	}
}; // InverseArrayIterator

static inline int nextDelta(int sz, int _delta)
{
	return (sz%_delta) ? ((sz+_delta)/_delta)*_delta : sz;
} // nextDelta

/*-------------------------------------------- */
/* class Array                                */
/* A dynamic array of pointers to objects      */
/*-------------------------------------------- */

/* ---- Array ---- */
template<class T>
Array<T>::Array(const Array<T>& v) :
	_capacity(v._capacity),
	_count(v._count),
	_delta(v._delta),
	_compare(v._compare)
{
	data = new T[v._capacity];
	assert(_capacity == 0 || (data != NULL && v.data != NULL) );
	if(v._count > 0)
		copyData(v.data, v._count);
} /* Array */

/* ---- operator = ---- */
template<class T>
Array<T>& Array<T>::operator = (const Array<T>& v)
{
	if (data != v.data) {
		if (data) delete [] data;
		if (v.data) {
			data = new T[v._capacity];
			assert(data != NULL);
			copyData(v.data, v._count);
		} else
			data = NULL;

		_capacity = v._capacity;
	}
	_count  = v._count;
	_delta  = v._delta;
	_compare= v._compare;
	return *this;
} /* operator = */

/* ---- copyData ---- */
template<class T>
void Array<T>::copyData(T* src, int n)
{
	if (IsPointer(data[0]))
		memcpy(data, src, sizeof(T)*n);
	else
		for (int i=0; i<n; i++)
			data[i] = src[i];	// involve the copy constructor
} /* copyData */

/* ---- zero ---- */
template<class T>
void Array<T>::zero(int idx)
{
	if (IsPointer(data[0]))
		memset(data+idx, 0, sizeof(T));
} /* zero */

/* ---- zero ---- */
template<class T>
void Array<T>::zero(int lwr, int upr)
{
	if (upr>_capacity) upr=_capacity;
	if (IsPointer(data[0]))
		memset(data+lwr, 0, sizeof(T)*(upr-lwr));
} /* zero */

/* ---- resize ---- */
template<class T>
bool Array<T>::resize(int newSz)
{
	if (newSz <= _capacity) return true;
	if (_delta == 0) return false;
	int sz = _capacity + nextDelta(newSz - _capacity, _delta);

	T* oldData = data;
	data = new T[sz];
	if (!data) {
		delete [] oldData;
		// throw an exception
		return false;
	}
	if (oldData) {
		copyData(oldData, _count);
		delete [] oldData;
	}

	_capacity = sz;
	zero(_count, _capacity);
	return true;
} /* resize */

/* ---- add ---- */
template<class T>
int Array<T>::add(const T &t)
{
	assert(_count<=_capacity);
	assert(_delta>0);
	if (_compare==NULL || !_count) {
		if (_count >= _capacity) {
			if (_delta<=0 || !resize(_count+1))
				return -1;
			_delta <<= 1;	// double delta on next iteration
		}
		data[_count] = t;
		return _count++;
	} else {
		if (_compare(t, data[0]) <= 0) {
			insert(0, t);
			return 0;
		}
		int upper = top();
		if (_compare(t, data[upper]) >= 0) {
			if (_count >= _capacity) {
				resize(_count+1);
				_delta <<= 1;	// double delta on next iteration
			}
			data[_count] = t;
			return _count++;
		}
		int lower = 0;
		while (lower < upper) {
			int middle = (lower+upper)/2;
			int cmp = _compare(t, data[middle]);
			if (cmp == 0) {
				insert(middle, t);
				assert(middle==0     || _compare(data[middle-1], data[middle]) <= 0);
				assert(middle==top() || _compare(data[middle], data[middle+1]) <= 0);
				return middle;
			}
			if (cmp > 0)
				lower = middle+1;
			else
				upper = middle;
		}
		/* add at low position */
		insert(lower, t);

		assert(lower==0     || _compare(data[lower-1], data[lower]) <= 0);
		assert(lower==top() || _compare(data[lower], data[lower+1]) <= 0);

		return lower;
	}
} /* add */

/* --- shiftby --- */
template<class T>
void Array<T>::shiftby(int idx, int len)
{
	if (_count+len-1 >= _capacity) {
		resize(_count+len);
		_delta <<= 1;	// double delta on next iteration
	}
	memmove(data+(idx+len), data+idx, sizeof(T)*(_count-idx));
	_count += len;
} /* shiftby */

/* ---- replace ---- */
template<class T>
void Array<T>::replace(int idx, const T &t)
{
	assert(_compare == NULL);
	data[idx] = t;
} /* replace */

/* ---- erase ---- */
template<class T>
void Array<T>::erase(const T &t)
{
	int pos = findPtr(t);
	if (pos>=0)
		detach(pos);
	//else
	//	return NULL;
} /* erase */

/* ---- erase ---- */
template<class T>
void Array<T>::erase(int idx)
{
	assert(idx>=0 && idx<_count);
//	T t = data[idx];
	_count--;
	memmove(data+idx, data+(idx+1), sizeof(T)*(_count-idx));
#ifdef _DEBUG
	zero(_count);
#endif
	//return t;
} /* erase */

/* ---- erase ---- */
template<class T>
void Array<T>::erase(int lwr, int upr)
{
	if (lwr>=_count) return;
	assert(lwr>=0);
	assert(upr>lwr);

	if (upr>_count) upr=_count;
	int n = upr-lwr;
	memmove(data+lwr, data+upr, sizeof(T)*(_count-upr));
	_count -= n;
#ifdef _DEBUG
	zero(_count, _count+n);
#endif
} /* erase */

/* ---- clear ---- */
template<class T>
void Array<T>::clear( )
{
	zero(0, _count);
	_count = 0;
} /* clear */

/* ---- find ---- */
template<class T>
int Array<T>::find(const T &t) const
{
	if (_compare == NULL) {
		if (_count != 0) {
			for (int idx = 0; idx < _count; idx++ )
				if (data[idx] == t)
					return idx;
		}
		return -1;
	} else {
		int lower = 0;
		int upper = top();
		if (_count != 0) {
			while (lower < upper) {
				int middle = (lower+upper)/2;
				int cmp = _compare(t, data[middle] );
				if (cmp == 0)
					return middle;
				if (cmp > 0)
					lower = middle+1;
				else
					upper = middle-1;
			}
		}
		if (lower == upper && _compare(data[lower], t) == 0)
			return lower;
		else
			return -1;
	}
} /* find */

/* ---- search ---- */
template<class T>
int Array<T>::search(const T &t) const
{
	assert(_compare);
	if (_count==0) return -1;
	int lower = 0;
	int upper = top();

	if (_compare(t, data[0])<0)     return -1;
	if (_compare(t, data[upper])>0) return top()+1;

	while (lower < upper) {
		int middle = (lower+upper)/2;
		if (middle==lower) return middle;
		int cmp = _compare(t, data[middle] );
		if (cmp > 0)
			lower = middle;
		else
		if (cmp < 0)
			upper = middle;
		else
			return middle;
	}
	return -1;
} /* search */

/* ---- findPtr ---- */
template<class T>
int Array<T>::findPtr(const T &t) const
{
	for (int idx=0; idx < _count; idx++ )
		if (data[idx] == t)
			return idx;
	return -1;
} /* findPtr */

/* ---- sort ---- */
template<class T>
void Array<T>::sort()
{
	assert(_compare != NULL);

	for (int i=0; i<_count; i++) {
		int change = 0;
		for (int j=top(); j>i; j--) {
			int jj = j-1;
			if (_compare(data[j], data[jj]) < 0) {
				/* swap values */
				T tmp = data[jj];
				data[jj]  = data[j];
				data[j]   = tmp;
				change = 1;
			}
		}
		if (!change)
			return;
	}
} /* sort */

/* ---- moveToEnd ---- */
template<class T>
void Array<T>::moveToEnd(size_t idx)
{
	if (idx==top()) return;
	T element = get(idx);
	detach(idx);
	push(element);
} /* moveToEnd */

/* ---- moveToEnd ---- */
template<class T>
void Array<T>::moveToEnd(T element)
{
	int idx = find(element);
	if (idx < 0) return;
	if (idx==top()) return;
	detach(idx);
	push(element);
} /* moveToEnd */

#endif
