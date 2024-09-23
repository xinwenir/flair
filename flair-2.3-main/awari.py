# 
#
# Copyright and User License
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# Copyright 2006-2019 CERN and INFN
# 
#
# Please consult the LICENSE file for the license 
#
# DISCLAIMER
# ~~~~~~~~~~
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, IMPLIED WARRANTIES OF MERCHANTABILITY, OF
# SATISFACTORY QUALITY, AND FITNESS FOR A PARTICULAR PURPOSE
# OR USE ARE DISCLAIMED. THE COPYRIGHT HOLDERS AND THE
# AUTHORS MAKE NO REPRESENTATION THAT THE SOFTWARE AND
# MODIFICATIONS THEREOF, WILL NOT INFRINGE ANY PATENT,
# COPYRIGHT, TRADE SECRET OR OTHER PROPRIETARY RIGHT.
#
# LIMITATION OF LIABILITY
# ~~~~~~~~~~~~~~~~~~~~~~~
# THE COPYRIGHT HOLDERS AND THE AUTHORS SHALL HAVE NO
# LIABILITY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
# CONSEQUENTIAL, EXEMPLARY, OR PUNITIVE DAMAGES OF ANY
# CHARACTER INCLUDING, WITHOUT LIMITATION, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES, LOSS OF USE, DATA OR PROFITS,
# OR BUSINESS INTERRUPTION, HOWEVER CAUSED AND ON ANY THEORY
# OF CONTRACT, WARRANTY, TORT (INCLUDING NEGLIGENCE), PRODUCT
# LIABILITY OR OTHERWISE, ARISING IN ANY WAY OUT OF THE USE OF
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGES.
#
# Author:	Vasilis.Vlachoudis@cern.ch
# Date:	13-Sep-2008


__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import sys
import math
import random
import string
try:
	from tkinter import *
except ImportError:
	from tkinter import *

import tkDialogs

_FILE = "awari.var"

_COLOR_BOARD    = "gray"
_COLOR_INACTIVE = "darkgray"
_COLOR_YOU      = "darkblue"
_COLOR_CPU      = "darkred"
_COLOR_SCORE    = "darkgreen"
_COLOR_ACTIVE   = "lightgreen"
_COLOR_SELECT   = "lightgreen"
_COLOR_MOVE     = "lightgray"


toplevel = None
awari    = None
youwin   = 0
cpuwin   = 0

#===============================================================================
# Awari Game
#===============================================================================
class Awari(Frame):
	"""
                                 A  W  A  R  I                        o o
                               ~~~~~~~~~~~~~~~~~                    ____oo_
                                                         1989-2009 /||    |\\
                                   My Side                 Vasilis  ||    |
                                                                    `.___.'
                           | 6 | 5 | 4 | 3 | 2 | 1 |
                           |ooo|ooo|ooo|ooo|ooo|ooo|
                  My home  |   |   |   |   |   |   | Your home
                           |ooo|ooo|ooo|ooo|ooo|ooo|
                           | 1 | 2 | 3 | 4 | 5 | 6 |

                                 Your  Side

             Awari is an  African  game  played  with seven sticks and
          thirty-six stones or  beans  laid  out  as  shown above. The
          board is divided into six compartments or pits on each side.
          In addition, there are two special HOME pits at the end.

             A  move  is  made  by  taking  all  the  beans  from  any
          (non-empty) pit on your own  side.  Starting from the pit to
          the right of this one,  these  beans  are "sown" one in each
          pit working around the board anticlockwise.

             A turn consists of one or  two moves. If the last bean of
          your move is sown in  your  own  home  you may take a second
          move.

             If the last bean is sown in a move lands in an empty pit,
          provided that the opposite pit is no empty, all the beans in
          the opposite  pit,  together  with  the  last  bean sown are
          "captured" and moved to the players home.

             When either side  is  empty,  the  game  is finished. The
          player with the most beans in his home has won.

             There is a learning mechanism  in the program that causes
          the play of the computer to improve as it plays more games.

                      email: Vasilis.Vlachoudis@cern.ch
	"""

	# ----------------------------------------------------------------------
	def __init__(self, master, imgDir, *arg, **kw_args):
		Frame.__init__(self, master, *arg, **kw_args)
		self.canvas = Canvas(self, background=_COLOR_BOARD)
		canvas = self.canvas
		canvas.pack(side=TOP)

		self.helpframe = Frame(self)
		self.helpframe.pack(side=BOTTOM, expand=YES, fill=BOTH)
		help = Text(self.helpframe)
		help.pack(side=LEFT, expand=YES, fill=BOTH)
		help.insert(0.0, self.__doc__)
		sb = Scrollbar(self.helpframe, orient=VERTICAL)
		sb.pack(side=RIGHT, fill=Y)
		sb.config(command=help.yview)
		help.config(yscrollcommand=sb.set)
		help.config(state=DISABLED)
		self.helpframe.forget()

		self.learn = []
		self.yourturn = True

		self.stick1 = PhotoImage(file = imgDir+"/stick1.gif")
		self.stick2 = PhotoImage(file = imgDir+"/stick2.gif")
		self.bean_image = []
		for i in range(1,5):
			self.bean_image.append(PhotoImage(
				file = "%s/red_bean%d.gif"%(imgDir,i)))

		height = self.stick1.height()

		self.margin     = 20		# margin
		self.beansize   = self.bean_image[0].width()
		self.space      = height * 3/8 + self.beansize/3	# dimension of pits
		self.pits       = [None]*14	# my pits
		self.bean       = []		# bean images
		self.board_text = [0]*14	# board text
		self.board_beans= [None]*14	# board beans links
		self.position   = [None]*14	# center of pits

		# fill board beans (double list)
		for i in range(14):
			self.board_beans[i] = []

		# --- Overall space ---
		canvas_width  = 8*self.space  + 2*self.margin

		# --- Text ---
		canvas.create_text(canvas_width/2+2, 2+2,
				text="AWARI",
				anchor=N, font=("Times",-40,"bold"),
				fill="darkgray")
		t = canvas.create_text(canvas_width/2, 2,
				text="AWARI",
				anchor=N, font=("Times",-40,"bold"),
				fill="darkred")
		x1,y1,x2,y2 = canvas.bbox(t)

		t2 = canvas.create_text(canvas_width/2, y2-3,
				text=__email__,
				anchor=N, font=("Times",-12,"italic"),
				fill="darkred")
		x1,y1,x2,y2 = canvas.bbox(t2)

		self.topmargin = y2+y2-y1+2		# top margin
		canvas_height = self.topmargin + 2*self.margin + height
		canvas.config(width=canvas_width, height=canvas_height)

		# ---
		font = ("Helvetica", -24)
		self.your_turn = canvas.create_text(
				self.margin + self.space, canvas_height,
				text="Your turn",
				anchor=SW, font=font,
				fill=_COLOR_YOU)

		self.cpu_turn = canvas.create_text(
				canvas_width - self.margin - self.space,
				canvas_height,
				text="My turn",
				anchor=SE, font=font,
				fill=_COLOR_INACTIVE)

		font = ("Helvetica", -18)
		self.message = canvas.create_text(
				canvas_width/2, canvas_height,
				text="You x - x Cpu",
				fill=_COLOR_SCORE,
				justify=CENTER,
				anchor=S, font=font)

		self.helptxt = canvas.create_text(canvas_width, canvas_height,
				text="Help",
				font=font,
				anchor=SE, fill="DarkOrange")
		canvas.tag_bind(self.helptxt, "<Button-1>", self.help)


		# ---------- sticks and pits -------
		x = self.margin
		y = height / 2 + self.topmargin
		xn = x + self.space

		self.pits[13] = canvas.create_rectangle((x,y-height/4,xn,y+height/4),
				fill=_COLOR_BOARD, width=0)

		font = ("Helvetica", -14)
		canvas.create_text(x+self.space/2, y+height/4-1, text="My home",
				fill=_COLOR_INACTIVE, font=font, anchor=S)

		y1 = y - height/2+3
		y2 = y - 3
		y3 = y + 3
		y4 = y + height/2-3

		xs = xn
		x  = xn
		xn = x + self.space
		self.pits[0] = canvas.create_rectangle((x,y3,xn,y4),
				fill=_COLOR_BOARD, width=0)
		self.pits[12] = canvas.create_rectangle((x,y1,xn,y2),
				fill=_COLOR_BOARD, width=0)
		canvas.create_image(x,y,image=self.stick1)

		x  = xn
		xn = x + self.space
		self.pits[1] = canvas.create_rectangle((x,y3,xn,y4),
				fill=_COLOR_BOARD, width=0)
		self.pits[11] = canvas.create_rectangle((x,y1,xn,y2),
				fill=_COLOR_BOARD, width=0)
		canvas.create_image(x,y,image=self.stick2)

		x  = xn
		xn = x + self.space
		self.pits[2] = canvas.create_rectangle((x,y3,xn,y4),
				fill=_COLOR_BOARD, width=0)
		self.pits[10] = canvas.create_rectangle((x,y1,xn,y2),
				fill=_COLOR_BOARD, width=0)
		canvas.create_image(x,y,image=self.stick1)

		x  = xn
		xn = x + self.space
		self.pits[3] = canvas.create_rectangle((x,y3,xn,y4),
				fill=_COLOR_BOARD, width=0)
		self.pits[9] = canvas.create_rectangle((x,y1,xn,y2),
				fill=_COLOR_BOARD, width=0)
		canvas.create_image(x,y,image=self.stick2)

		x  = xn
		xn = x + self.space
		self.pits[4] = canvas.create_rectangle((x,y3,xn,y4),
				fill=_COLOR_BOARD, width=0)
		self.pits[8] = canvas.create_rectangle((x,y1,xn,y2),
				fill=_COLOR_BOARD, width=0)
		canvas.create_image(x,y,image=self.stick1)

		x  = xn
		xn = x + self.space
		self.pits[5] = canvas.create_rectangle((x,y3,xn,y4),
				fill=_COLOR_BOARD, width=0)
		self.pits[7] = canvas.create_rectangle((x,y1,xn,y2),
				fill=_COLOR_BOARD, width=0)
		canvas.create_image(x,y,image=self.stick2)

		x  = xn
		xn = x + self.space
		canvas.create_image(x,y,image=self.stick1)

		self.pits[6] = canvas.create_rectangle((x,y-height/4,xn,y+height/4),
				fill=_COLOR_BOARD, width=0)
		canvas.create_text(x+self.space/2, y+height/4-1, text="Your home",
				fill=_COLOR_INACTIVE, font=font, anchor=S)

		canvas.tag_lower(canvas.create_line(xs,y, x,y, fill="darkred"))
		for r in self.pits:
			if r: canvas.addtag_withtag("pits",r)
		canvas.tag_lower("pits")

		width = x + self.space

		# ---------- numbers ---------
		x = self.margin + self.space/2
		y = self.topmargin + height
		font = ("Helvetica", -14)
		for i in range(6):
			x += self.space
			self.board_text[i] = canvas.create_text(x,y,text=str(i+1),
					anchor=N, font=font, fill=_COLOR_YOU)
			self.position[i] = (x, y-self.space/2)

			self.board_text[12-i] = canvas.create_text(x,self.topmargin,
					text=str(6-i),
					anchor=S, font=font ,fill=_COLOR_CPU)
			self.position[12-i] = (x, self.topmargin+self.space/2)

		y = self.topmargin + height/2

		self.board_text[ 6] = canvas.create_text(canvas_width-self.margin/2, y,
				text="m",anchor=E, font=font, fill=_COLOR_YOU)
		self.position[ 6] = (width-self.space/2,y)

		self.board_text[13] = canvas.create_text(self.margin/2, y,
				text="y",anchor=W, font=font, fill=_COLOR_CPU)
		self.position[13] = (self.margin+self.space/2,y)

		# -------- create beans --------
		for i in range(6):
			for j in range(3):
				x,y = self.rndPosition(i)
				bi = random.randint(0,len(self.bean_image)-1)
				b = canvas.create_image(x,y,image=self.bean_image[bi])
				self.board_beans[i].append(b)
				self.bean.append(b)

				x,y = self.rndPosition(12-i)
				bi = random.randint(0,len(self.bean_image)-1)
				b = canvas.create_image(x,y,image=self.bean_image[bi])
				self.board_beans[12-i].append(b)
				self.bean.append(b)

		# ---- bring numbers up ---
		for t in self.board_text:
			self.canvas.tag_raise(t)

		self.canvas.bind("<Button-1>", self.click)

	# ----------------------------------------------------------------------
	def reset(self):
		self.yourturn = True
		self.end      = False
		self.again    = False
		self.board    = [3]*14; self.board[6] = 0; self.board[13] = 0
		self.move     = 0

		for i in range(14):
			self.board_beans[i] = []
		k = 0
		for i in range(6):
			for j in range(3):
				self.moveBean(self.bean[k], i)
				k += 1
				self.moveBean(self.bean[k], 12-i)
				k += 1
		self.updateBoard()
		self.canvas.itemconfig("pits", fill=_COLOR_BOARD)
		self.learn.append(0)
		self.activatePits()

	# ----------------------------------------------------------------------
	def read(self, filename):
		global cpuwin, youwin
		del self.learn[:]
		self.filename = filename
		try:
			awarif = open(self.filename, "r")
			youwin, cpuwin = list(map(int,awarif.readline().split()))
			self.board = list(map(int,awarif.readline().split()))
			self.learn = list(map(int, awarif.readlines()))
			awarif.close()
		except:
			pass

	# ----------------------------------------------------------------------
	def write(self):
		# Save everything we done till now
		try:
			awarif = open(self.filename, "w")
			awarif.write("%d %d\n"%(youwin, cpuwin))
			awarif.write("%s\n"%(" ".join(map(str,self.board))))
			for i in self.learn:
				if i!=0: awarif.write("%d\n"%(i))
			awarif.close()
		except:
			pass

	# ----------------------------------------------------------------------
	# Find best move for CPU
	# ----------------------------------------------------------------------
	def CPU_Move(self):
		D = -99
		Q = -99
		G = self.board[:]		# Save the board
		h = 13

		# Calculate all possible moves
		#self.printBoard()
		for j in range(7,13):
			if self.board[j] == 0: continue
			self.moveIt(j, h, False)

			for i in range(6):
				if self.board[i] == 0: continue
				L = self.board[i] + i	# last
				R = (L>13)		# full turn
				L = L % 14
				P = 12 - L		# opposite
				if self.board[L]==0 and L!=6 and L!=13:	# Again the same story
					R += self.board[P]
				if R > Q: Q = R

			Q = self.board[13] - self.board[6] - Q
			if self.move < 9:
				k = j
				if k > 6: k -= 7
				for lea in self.learn[:-1]:
					if self.learn[-1]*6+k == int(lea/pow(6,
								(7-self.move)+0.1)):
						Q -= 2

			for i in range(14):
				self.board[i] = G[i]	# Restore the board

			if Q >= D:	# if the current move is better
				A = j	# then save it
				D = Q

		m = A
		return A-6, self.makeMove(m, h)

	# ----------------------------------------------------------------------
	# Perform the move
	# ----------------------------------------------------------------------
	def makeMove(self, m, h):
		k = m
		m = self.moveIt(m, h, True)
		self.end = False
		if k > 6: k -= 7
		self.move += 1

		# self.learn dummit
		if self.move < 9:
			self.learn[-1] = self.learn[-1]*6 + k

		# Game Over?
		for i in range(7,13):
			if self.board[i] != 0:
				for j in range(6):
					if self.board[j] != 0:
						return m
		self.end = True
		return m	# return last pit

	# ----------------------------------------------------------------------
	# Move the 'M' pit
	# ----------------------------------------------------------------------
	def moveIt(self, m, h, update):
		if update:
			pitbeans = self.board_beans[m]
			self.board_beans[m] = []
		p = self.board[m]
		self.board[m] = 0
		for i in range(p):
			m += 1
			if m > 13: m -= 14
			self.board[m] += 1
			if update:
				self.moveBean(pitbeans[i], m)

		# if last pit was empty and not home
		p = 12 - m
		if self.board[m]==1 and m!=6 and m!=13:
			# and opposite pit was empty too
			if self.board[p] != 0:
				if update:
					b = self.board_beans[m].pop()
					self.moveBean(b, h)
					for b in self.board_beans[p]:
						self.moveBean(b, h)
					self.board_beans[p] = []
					self.canvas.itemconfig(self.pits[p], fill=_COLOR_MOVE)
				#
				self.board[h] += self.board[p] + 1
				self.board[m] = 0
				self.board[p] = 0
		if update:
			self.updateBoard()
		return m	# return last pit

	# ----------------------------------------------------------------------
	# Procedure for printing the board
	# ----------------------------------------------------------------------
	def printBoard(self):
		sys.stdout.write()
		cc = ' '*22
		p = cc+'  '
		for i in range(12,6,-1):
			p += "%4d"%(self.board[i])
		sys.stdout.write("%s\n"%(p))
		sys.stdout.write("%s %2d                          %d\n" % \
				(cc,self.board[13],self.board[6]))
		p = cc+'  '
		for i in range(6):
			p += "%4d"%(self.board[i])
		sys.stdout.write("%s\n"%(p))

	# ----------------------------------------------------------------------
	def moveBean(self, bean, i):
		x,y = self.rndPosition(i)
		self.canvas.coords(bean, x, y)
		self.board_beans[i].append(bean)
		self.canvas.itemconfig(self.pits[i], fill=_COLOR_MOVE)

	# ----------------------------------------------------------------------
	def updateBoard(self):
		for j in range(14):
			self.canvas.itemconfig(self.board_text[j],
				text=str(len(self.board_beans[j])))
		self.canvas.itemconfig(self.message,
				text="You %d - %d Cpu"%(youwin,cpuwin))

	# ----------------------------------------------------------------------
	def activatePits(self):
		if self.yourturn:
			for i in range(6):
				if self.board[i] > 0:
					self.canvas.itemconfig(self.pits[i],
						activefill=_COLOR_ACTIVE)
		else:
			self.canvas.itemconfig("pits", activefill="")

	# ----------------------------------------------------------------------
	def rndPosition(self, idx):
		for i in range(20):
			x,y = self.position[idx]
			r   = random.randint(0, self.space/3)
			phi = random.random() * 2.0 * math.pi
			x += r * math.cos(phi)
			y += r * math.sin(phi)

			# compare for overlap against other beans
			ok = True
			for b in self.board_beans[idx]:
				xb, yb = self.canvas.coords(b)
				if (x-xb)**2 + (y-yb)**2 < self.beansize**2:
					ok = False
					continue
			if ok: break
		return (x,y)

	# ----------------------------------------------------------------------
	def click(self, event):
		if not self.yourturn: return
		x = int(self.canvas.canvasx(event.x))
		y = int(self.canvas.canvasy(event.y))
		for i in self.canvas.find_overlapping(x,y,x+1,y+1):
			try:
				m = self.pits.index(i)
			except:
				break
			if m<6 and self.board[m]>0:
				self.canvas.itemconfig("pits", fill=_COLOR_BOARD)
				self.canvas.itemconfig(i, fill=_COLOR_SELECT)
				if self.again:
					txt = self.canvas.itemcget(self.your_turn,"text")
					self.canvas.itemconfig(self.your_turn,
						text="%s, %d"%(txt,m+1))
				else:
					self.canvas.itemconfig(self.your_turn,
						text="Your turn: %d"%(m+1))
				m = self.makeMove(m, 6)
				if self.end:
					self.yourturn = False
					self.activatePits()
					self.gameOver()
					return
				#self.printBoard()
				if m==6 and not self.again:
					self.canvas.itemconfig(self.message,
						text="play again")
					self.again = True
				else:
					self.again = False
					self.yourturn = False
					self.activatePits()
					self.canvas.itemconfig("pits", activefill="")
					self.canvas.itemconfig(self.your_turn,fill=_COLOR_INACTIVE)
					self.canvas.itemconfig(self.cpu_turn,fill=_COLOR_CPU,
								text="CPU turn:")
					self.after(2000, self._timer)

	# ----------------------------------------------------------------------
	def _timer(self):
		self.canvas.itemconfig("pits", fill=_COLOR_BOARD)
		if self.end:
			self.gameOver()
			return
		m1,l = self.CPU_Move()
		self.canvas.itemconfig(self.cpu_turn,
				text="CPU turn: %d"%(m1))
		self.canvas.itemconfig(self.pits[m1+6], fill=_COLOR_SELECT)
		if self.end:
			self.gameOver()
			return
		if l == 13:
			m2,l = self.CPU_Move()
			self.canvas.itemconfig(self.pits[m2+6], fill=_COLOR_SELECT)
			self.canvas.itemconfig(self.cpu_turn,
					text="CPU turn: %d, %d"%(m1,m2))
			if self.end:
				self.gameOver()
				return
		#self.printBoard()

		self.yourturn = True
		self.activatePits()
		self.canvas.itemconfig(self.your_turn, fill=_COLOR_YOU,
					text="Your turn:")
		self.canvas.itemconfig(self.cpu_turn,fill=_COLOR_INACTIVE)

	# ----------------------------------------------------------------------
	def gameOver(self):
		global cpuwin, youwin
		msg = "Game Over\n"
		d = self.board[6] - self.board[13]
		if d < 0:
			Sa =  "I win by %d point"%(-d)
			cpuwin += 1
			self.learn.pop()
		elif d == 0:
			Sa = "Drawn game"
		else:
			Sa =  "You win by %d point"%(d)
			youwin += 1

		if abs(d)<=1:
			msg += "%s."%(Sa)
		else:
			msg += "%ss."%(Sa)

		self.updateBoard()
		self.canvas.itemconfig(self.message, text=msg)

		ans = tkDialogs.Dialog(toplevel,
			{"title":"Awari: Game over",
			 "text":msg,
			 "bitmap":"questhead",
			 "default":0,
			 "strings":('New game', 'Close')})

		if ans.num == 0:
			self.reset()
		else:
			close()

	# ----------------------------------------------------------------------
	def help(self, event):
		if self.helpframe.winfo_ismapped():
			self.canvas.itemconfig(self.helptxt, text="Help")
			self.helpframe.forget()
		else:
			self.canvas.itemconfig(self.helptxt, text="Hide Help")
			self.helpframe.pack(side=BOTTOM)

#-------------------------------------------------------------------------------
def close():
	global toplevel, awari
	if not awari.end:
		awari.learn.pop()
	toplevel.destroy()
	awari.write()
	toplevel = None
	awari = None

#-------------------------------------------------------------------------------
def startGame(master, filename, imgDir):
	global toplevel, awari
	if toplevel is None:
		toplevel = Toplevel(master)
		toplevel.title("awari")
		toplevel.protocol("WM_DELETE_WINDOW", close)
		awari = Awari(toplevel, imgDir)
		awari.pack(expand=YES, fill=BOTH)
		awari.read(filename)
		awari.reset()
	else:
		toplevel.deiconify()

#-------------------------------------------------------------------------------
if __name__ == "__main__":
	root = Tk()
	startGame(None,"awari.var", "icons")
	root.mainloop()
