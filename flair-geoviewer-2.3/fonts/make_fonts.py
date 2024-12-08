#!/usr/bin/python2
# -*- coding: latin1 -*-
# $Id$
#
# Author:	Vasilis.Vlachoudis@cern.ch
# Date:	10-Dec-2012

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import os
import sys
import time
import random
import commands
#from Tkinter import *
import tkFont
import tkinter as tk
import tkinter.font as tkFont
def draw_canvas(font):
	width, height = metrics(font)
	width  += 1
	height += 1
	canvas.delete(ALL)
	print ("Font:", width,height)
	print ("Size:", width*16, height*16)
	tk.geometry("%dx%d"%(width*16+2,height*16+2))
	tk.update_idletasks()
	tk.lift()
	for j in range(16):
		y = j*height+1
		for i in range(16):
			x = i*width+1
			canvas.create_text(x,y,
				text=chr(j<<4 | i),
				font=font,
				fill="White",
				anchor=NW)
	tk.lift()
	tk.update_idletasks()
	time.sleep(1)

def metrics(font):
	width = 0
	for i in range(256):
		width = max(width, font.measure(chr(i)))
	return width, font.metrics()["linespace"]

def capture(name):
	print ("Canvas:", canvas.winfo_width(), canvas.winfo_height())
	window = None
	for line in commands.getoutput("wmctrl -l").splitlines():
		words = line.split()
		try:
			if words[3] == title:
				window = words[0]
				break
		except IndexError:
			pass

	if window is None: return

	os.system("xwd -out font.xwd -id %s -nobdrs"%(window))
	os.system("convert font.xwd -crop %dx%d+1+1 -orient TopLeft -compress none -flip -type Grayscale %s"%(canvas.winfo_width()-2, canvas.winfo_height()-2, name))

	f = open(name,"rb+")
	f.seek(17)
	f.write(chr(0x20))
	f.close()

def drawFont(family, size):
	print
	print ("Font:",family, size)
	font = tkFont.Font(tk,(family, size))
	draw_canvas(font)
	capture("%s-%02d.tga"%(family.replace(" ","_"),size))

def process(event=None):
	for family in ["Arial", "Helvetica", "Courier", "DejaVu Sans", "DejaVu Serif"]:
		for size in [8,9,10,11,12,14,16,18,20,22,24,28,32,36]:
			drawFont(family,size)

def quit(event=None):
	tk.destroy()

tk = Tk()
title = "Font.Capture-%.10g"%(random.random())
tk.title(title)
tk.bind("<Escape>", quit)
tk.geometry("+0+0")
tk.lift()

canvas = Canvas(tk, background="Black")
canvas.pack(fill=BOTH, expand=YES, padx=0, pady=0)
tk.wait_visibility()

tk.after(100,process)
tk.mainloop()
