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

#ifndef __LIST_H
#define __LIST_H

#include <assert.h>
#include <stdlib.h>

template <class T> class List;
template <class T> class ListIterator;

template <class T>
class ListElement {
private:
	T		data;
	ListElement<T>*	prev;
	ListElement<T>*	next;
public:
	ListElement() : prev(NULL), next(NULL) {}
	ListElement(T& d) : data(d), prev(NULL), next(NULL) {}
friend	class List<T>;
friend	class ListIterator<T>;
}; // class ListElement<T>

template <class T>
class List {
private:
	int		items;
	ListElement<T>*	head;
	ListElement<T>*	tail;
public:
	List(): items(0), head(NULL), tail(NULL) {}
	~List() {
			ListElement<T> *cur = head;
			while (cur) {
				ListElement<T> *nxt = cur->next;
				delete cur;
				cur = nxt;
			}
		}

	/** @return number of elements in the list */
	int	count() const { return items; }

	/** add an element at the end of the list
	 * @param dat data pointer to insert
	 */
	void	append(T dat) {
			items++;
			ListElement<T>* elem = new ListElement<T>(dat);
			elem->prev = tail;
			if (tail) tail->next = elem;
			tail = elem;
			if (head == NULL) head = elem;
		}

	List<T>& operator <<(T dat)	{ append(dat); return *this; }

	/** add an element as a head of the list
	 * @param dat data pointer to insert
	 */
	void	addHead(T dat) {
			items++;
			ListElement<T> *elem = new ListElement<T>(dat);
			elem->next = head;
			if (head) head->prev = elem;
			head = elem;
			if (tail == NULL) tail = elem;
		}

	/** append all elements from another list
	 * @param list list to append
	 */
	void	extend(List<T>* list) {
			ListElement<T>* cur = list->tail;
			while (cur) {
				append(cur->data);
				cur = cur->next;
			}
		}

	/** insert an element after the specified position
	 * @param index of position where to insert the element
	 * @param dat data pointer to insert
	 */
	void	insert(int index, T dat) {
			if (index == 0)
				addHead(dat);
			else
			if (index >= items)
				append(dat);
			else {
				ListElement<T> *prev = element(index);
				ListElement<T> *next = prev->next;
				ListElement<T> *elem = new ListElement<T>(dat);
				items++;
				elem->prev = prev;
				prev->next = elem;
				elem->next = next;
				next->prev = elem;
			}
		}

	/** detach an element from the list given by its index
	 * @param index of element to detach
	 */
	void	detach(int index) {
			ListElement<T> *elem = element(index);
			if (elem) delElement(elem);
		}

	/** detach element from the list by its pointer
	 * @param dat pointer of element to detach */
	void	detach(T dat) {
			ListElement<T> *elem = element(dat);
			if (elem) delElement(elem);
		}

//	void	reverse();
//	void	sort();

	/** find element by pointer
	 * @param dat element to find
	 * @return index number in list, -1 for non existent */
	int	find(T dat) {
			int idx = 0;
			ListElement<T> *cur = head;
			while (cur) {
				if (cur->data == dat) return idx;
				cur = cur->next;
				idx++;
			}
			return -1;
		}

	/** return and remove last element of the list
	 * @return last element of the list
	 */
	T	pop() {
			if (tail) {
				T dat = tail->dat;
				items--;
				ListElement<T> *prev = tail->prev;
				if (prev) prev->next = NULL;
				delete tail;
				tail = prev;
				if (tail==NULL) head=NULL;
				return dat;
			} else
				return NULL;
		}

	/** return an element in list given by its index
	 * @param index of element to be returned
	 * @return data pointer
	 */
	T	get(int index) const {
			ListElement<T>* elem = element(index);
			if (elem == NULL) return NULL;
			return elem->dat;
		}

	/**
	 * @param index index of data to be returned
	 * @return data pointer */
	T	operator [] (int index)	const	{ return get(index); }

	/** @return first element in list */
	T	first()				{ return head?head->dat:NULL; }

	/** @return last element in list */
	T	last()				{ return tail?tail->dat:NULL; }

private:
	ListElement<T>*	element(int index) const {
			ListElement<T> *cur = head;
			while (cur) {
				if (index==0) return cur;
				cur = cur->next;
				index--;
			}
			assert(index==0);
			return NULL;
		}

	ListElement<T>*	element(T dat) const {
			ListElement<T> *cur = head;
			while (cur) {
				if (cur->data == dat) return cur;
				cur = cur->next;
			}
			return NULL;
		}

	void	delElement(ListElement<T> *elem) {
			assert(elem != NULL);
			ListElement<T>* prev = elem->prev;
			ListElement<T>* next = elem->next;
			items--;
			delete elem;
			if (prev)
				prev->next = next;
			else
				head = next;
			if (next)
				next->prev = prev;
			else
				tail = prev;
		}

friend	class ListIterator<T>;
}; // class List<T>

/*------------------------------------------------------------------------*/
/** Iterator class for the List class
 * @see List
 */
template <class T>
class ListIterator {
private:
const	List<T>*	list;
	ListElement<T>*	cur;

public:
	/** Iterator constructor
	 * @param l List to scan
	 * @see List
	 */
	ListIterator(const List<T>& l) { reset(l); }

	/** reset list */
	void	reset(const List<T>& l) {
			list = &l;
			restart();
		}

	/** Test if end has reached
	 * @return true if the end has reached
	 */
	operator int()		const	{ return cur != NULL; }

	/** @return current data pointer */
	void*	current()		{ return cur->data; }

	/** Increase by one the iterator
	 * @return new data pointer
	 */
	T	operator ++(int) {
			T dat = cur->data;
			cur = cur->next;
			return dat;
		}

	/** Increase by one the iterator
	 * @return new data pointer
	 */
	T	operator ++() {
			void *dat = cur->data;
			cur = cur->next;
			return dat;
		}

	/** restart scanning */
	void	restart()	{ cur = list->head; }
}; // class ListIterator

#endif
