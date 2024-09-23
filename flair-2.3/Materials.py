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
# Date:	07-Mar-2007

__author__ = "Vasilis Vlachoudis"
__email__  = "Paola.Sala@mi.infn.it"

import sys
import time
import string

try:
	from tkinter import *
	import tkinter.messagebox as messagebox
except ImportError:
	from tkinter import *
	import tkinter.messagebox as messagebox

try:
	import configparser
except ImportError:
	import configparser as ConfigParser

import tkinter as tk
import undo
import Utils
import tkFlair
import tkExtra
import tkTextEditor
import PeriodicTable

import Input
import Ribbon
import Project
import FlairRibbon

#-------------------------------------------------------------------------------
# Global variables
#-------------------------------------------------------------------------------
groups    = []
materials = {}
changed   = False

_COMPOSITION = ["atom", "mass", "volume"]
_STATE       = ["solid", "liquid", "gas"]
_PAGE_TABS   = ["Notes", "Stoichiometry", "Properties"]

_clipboard = []

#-------------------------------------------------------------------------------
# Create material and insert to input
#-------------------------------------------------------------------------------
def insertMaterial(flair, mat):
	# Create the basic materials
	material_list = flair.project.input.materialList()

	# Create the material card and then the compound if needed
	matnames = mat.name.split()

	undoinfo  = [flair.refreshUndo(), flair.renumberUndo()]

	# Find a name that doesn't exist or replace material?
	for n in matnames:
		matname = n[:8]
		if matname=="PHOSPHOR": matname = "PHOSPHO"
		try:
			matid = material_list.index(matname)+1
		except ValueError:
			break
	else:
		#say("* Material",matnames,"already exists in input")
		matname = matnames[0][:8]
		if matname=="PHOSPHOR": matname = "PHOSPHO"
		if not tkFlair.askyesno("Material exists",
				"Material %s already exists!\n"\
				"Do you want to insert it anyway?"%(matname),
				parent=flair):
			return None, None

	# Add MATERIAL/COMPOUND cards after the last MATERIAL
	# or COMPOUND card
	pos = 0
	for tag in ["MATERIAL", "COMPOUND"]:
		try:
			card = flair.project.input.cards[tag][-1]
			pos = max(pos, card._pos)
		except:
			pass
	# If everything fails try GEOEND
	if pos == 0:
		try:
			card = flair.project.input.cards["GEOEND"][-1]
			pos = card._pos
		except:
			pass
	pos += 1

	tag = "MATERIAL"

	# Compound, or a single element?
	matelem = []
	if len(mat.elements)==1:
		# Single element
		Z, A, f = mat.elements[0]
		element = PeriodicTable.element(Z)
		if A == 0: A = ""

		what = [matname, Z, "", mat.density, "", "", A]
		card = Input.Card(tag, what)
		undoinfo.append(flair.addCardUndo(card, pos, False))
		pos += 1
#		say(card)
		matelem.append(matname)
	else:
		# Compound
		for Z, A, f in mat.elements:
			element = PeriodicTable.element(Z)
			if A==0:
				name  = element.name.upper()[:8]
				# Special case of Phosphorus
				if name=="PHOSPHOR": name = "PHOSPHO"
			else:
				name  = "%s-%d"%(element.symbol.upper(),A)
			if A == 0: A = ""

			matelem.append(name)
			try:
				# Check if the material exist
				material_list.index(name)
				#say("* Material",name,"already exists")
				continue
			except ValueError:
				pass

			what = [name, Z, "", element.density(), "", "", A]
			card = Input.Card(tag, what)
			undoinfo.append(flair.addCardUndo(card, pos, False))
			pos += 1

			# Check if the material exists?
#			say(card)

		what = [matname, "", "", mat.density]
		comment = mat.title + "\n" + mat.notes
		card = Input.Card(tag, what, comment)
		undoinfo.append(flair.addCardUndo(card, pos, False))
		pos += 1
#		say(card)

		# Create the compound
		tag = "COMPOUND"
		if   mat.composition == 0:
			s1 = 1.0
			s2 = ""
		elif mat.composition == 1:
			s1 = -1.0
			s2 = ""
		else:
			s1 = -1.0
			s2 = ""

#			for i in range(0,len(mat.elements),3):
#				what = [matname]
#				for j in range(3):
#					try:
#						Z, A, f = mat.elements[i+j]
#					except:
#						break
#					element = PeriodicTable.element(Z)
#					if A == 0:
#						ename = element.name.upper()[:8]
#					else:
#						ename = "%s-%d"%(element.symbol.upper(),A)
#					what.append(s1*f)
#					what.append(s2+ename)
#				card = Input.Card(tag, what)
#				undoinfo.append(flair.addCardUndo(card, pos, False))
#				pos += 1

		what = [matname]
		for i in range(0,len(mat.elements)):
			try: Z, A, f = mat.elements[i]
			except: break
			element = PeriodicTable.element(Z)
			if A == 0:
				ename = element.name.upper()[:8]
			else:
				ename = "%s-%d"%(element.symbol.upper(),A)
			ename = ename[:8]
			if ename=="PHOSPHOR": ename = "PHOSPHO"
			what.append(s1*f)
			what.append(s2+ename)
		remain = len(mat.elements)%3
		if remain > 0:
			while remain < 3:
				what.append("")
				what.append("")
				remain += 1
		card = Input.Card(tag, what)
		undoinfo.append(flair.addCardUndo(card, pos, False))
		pos += 1

	# Check for special material properties
	prop = []

	# Check for MAT-PROP
	if mat.pressure != "" or mat.rhor != "" or mat.ionization != "":
		prop.append("MAT-PROP")
	if mat.dpa_threshold != "":
		prop.append("DPA")
	if mat.s_cbar != "" or mat.s_x0 != "" or mat.s_x1 != "" or \
	   mat.s_a != "" or mat.s_m != "" or mat.s_d0 != "":
		prop.append("STERNHEI")
	if mat.xrho_dEdx != "" or mat.xrho != "":
		prop.append("CORRFACT")

	if prop:
		dlg = SelectPropertiesDlg(flair, matname, prop)
		if dlg.add:
			if "MAT-PROP" in dlg.add:
				what = ["", mat.pressure, mat.rhor,
					mat.ionization, matname]
				card = Input.Card("MAT-PROP", what)
				undoinfo.append(flair.addCardUndo(card, pos, False))
				pos += 1

			if "DPA" in dlg.add:
				what = ["DPA-ENER", mat.dpa_threshold,
					"", "", matname]
				card = Input.Card("MAT-PROP", what)
				undoinfo.append(flair.addCardUndo(card, pos, False))
				pos += 1

			if "STERNHEI" in dlg.add:
				what = [matname, mat.s_cbar, mat.s_x0,
					mat.s_x1, mat.s_a, mat.s_m,
					mat.s_d0]
				card = Input.Card("STERNHEI", what)
				undoinfo.append(flair.addCardUndo(card, pos, False))
				pos += 1

			if "CORRFACT" in dlg.add:
				what = ["", mat.xrho_dEdx, mat.xrho,
					"", matname]
				card = Input.Card("CORRFACT", what)
				undoinfo.append(flair.addCardUndo(card, pos, False))
				pos += 1

	undoinfo.append(flair.renumberUndo())
	undoinfo.append(flair.refreshUndo())

	flair.addUndo(undo.createListUndo(undoinfo, "Materials added"))

	# Return
	return matname, matelem

#-------------------------------------------------------------------------------
# Insert the material list to input
#-------------------------------------------------------------------------------
def insert2Input(flair, mats):
	if isinstance(mats,list):
		matnames = []
		elements = set()
		for m in mats:
			n, elems = insertMaterial(flair, materials[m])
			if n:
				matnames.append(n)
				for e in elems:
					elements.add(e)
		matnames = " ".join(matnames)
		elements = " ".join(sorted(list(elements)))
	else:
		matnames,elements = insertMaterial(flair, materials[mats])
		if elements: elements = " ".join(sorted(elements))

	if matnames:
		flair.notify("Material: %s"%(matnames),
			"Added materials:\n%s\n\n"
			"Based on:\n%s" % (matnames, elements),
			tkFlair.NOTIFY_NORMAL)
	else:
		flair.notify("No material added",
			"No material added to input.",
			tkFlair.NOTIFY_WARNING)
	return matnames, elements

#===============================================================================
# Material class
#===============================================================================
class Material:
	def __init__(self, title):
		self.title         = title
		self.name          = ""
		self.group         = []
		self.density       = 0.0
		self.composition   = 0
		self.notes         = ""
		self.elements      = []
		self.stoichiometry = ""

		self.state         = 0
		self.pressure      = ""	# Pressure
		self.rhor          = "" # Density factor
		self.ionization    = "" # Average Ionization potential
		self.dpa_threshold = "" # DPA Damage threshold
		self.xrho_dEdx     = ""
		self.xrho          = ""
		self.s_cbar        = ""	# Sternheimer parameters
		self.s_x0          = ""
		self.s_x1          = ""
		self.s_a           = ""
		self.s_m           = ""
		self.s_d0          = ""

		self.modified      = ""

	# ----------------------------------------------------------------------
	def make(self):
		# Make stoichiometry
		self.makeStoichiometry()

		# Calculate chksum
		self.chksum  = 0.0
		try: self.chksum  = self.density
		except: pass

		try: self.chksum += float(self.pressure)
		except: pass
		try: self.chksum += float(self.rhor)
		except: pass
		try: self.chksum += float(self.ionization)
		except: pass
		try: self.chksum += float(self.dpa_threshold)
		except: pass
		try: self.chksum += float(self.xrho_dEdx)
		except: pass
		try: self.chksum += float(self.xrho)
		except: pass

		try: self.chksum += float(self.s_cbar)
		except: pass
		try: self.chksum += float(self.s_x0)
		except: pass
		try: self.chksum += float(self.s_x1)
		except: pass
		try: self.chksum += float(self.s_a)
		except: pass
		try: self.chksum += float(self.s_m)
		except: pass
		try: self.chksum += float(self.s_d0)
		except: pass

		for Z,A,f in self.elements:
			try: self.chksum += float(Z)
			except: pass
			try: self.chksum += float(A)
			except: pass
			try: self.chksum += float(f)
			except: pass

	# ----------------------------------------------------------------------
	def setModified(self):
		self.modified = time.strftime("%d/%m/%Y",time.localtime())

	# ----------------------------------------------------------------------
	# Set the stoichiometry string
	# ----------------------------------------------------------------------
	def makeStoichiometry(self):
		if len(self.elements)==1:
			Z, A, f = self.elements[0]
			if A==0: A=""
			else: A=str(A)
			element = PeriodicTable.element(Z)
			self.stoichiometry = " %s%s" % (A, element.symbol)
			return

		# Percentage or absolute
		mult = 1.0
		try:
			if max([f for Z, A, f in self.elements])<1.0:
				mult = 100.0
		except ValueError:
			pass

		cnt = 0
		stoichiometry = ""
		for Z, A, f in self.elements:
			if A==0: A=""
			else: A=str(A)

			f *= mult
			if float(int(f)) == f:
				f = str(int(f))
			else:
				if mult>1.0:
					f = ("%.3g"%(f)).strip()
				else:
					f = str(f)

			element = PeriodicTable.element(Z)
			stoichiometry += " %s%s-%s," % (A, element.symbol, f)
			cnt += 1
			if cnt == 10: break

		if len(stoichiometry)>0:
			self.stoichiometry = stoichiometry[1:-1]
		else:
			self.stoichiometry = ""


#===============================================================================
# Material Frame
#===============================================================================
class MaterialPage(FlairRibbon.FlairPage):
	"""Database of predefined and user materials."""

	_name_ = "Materials"
	_icon_ = "materialdb"

	#----------------------------------------------------------------------
	def createRibbon(self):
		FlairRibbon.FlairPage.createRibbon(self)

		# ========== Input ===========
		group = Ribbon.LabelGroup(self.ribbon, "Materials")
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)
		group.grid3rows()

		# ---
		col,row = 0,0
		b = Ribbon.LabelButton(group.frame,
				image=tkFlair.icons["input32"],
				text="To Input",
				compound=TOP,
				command=self.insert2Input,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=3, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Insert selected materials to input")
		self.list_widgets.append(b)

		# ---
		col,row = 1,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["import32"],
				text="Import",
				compound=TOP,
				command=self.importFromInput,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, rowspan=3, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Import Materials FROM Input")
		self.list_widgets.append(b)

		# ---
		col,row = 2,0
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["add"],
				text="Add",
				compound=LEFT,
				anchor=W,
				command=self.add,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Create new material")
		self.list_widgets.append(b)

		# ---
		col,row = 2,1
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["x"],
				text="Delete",
				compound=LEFT,
				anchor=W,
				command=self.delete,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Clone run [Ctrl-D]")
		self.list_widgets.append(b)

		# ---
		col,row = 2,2
		b = Ribbon.LabelButton(group.frame, image=tkFlair.icons["edit"],
				text="Edit",
				compound=LEFT,
				anchor=W,
				command=self.edit,
				background=Ribbon._BACKGROUND)
		b.grid(row=row, column=col, padx=0, pady=0, sticky=NSEW)
		tkExtra.Balloon.set(b, "Edit material")
		self.list_widgets.append(b)

		# -------------- Edit Ribbon -------------
		self.editribbon = group = Ribbon.LabelGroup(self.ribbon, "Edit")
		group.label["background"] = Ribbon._BACKGROUND_GROUP2
		group.pack(side=LEFT, fill=Y, padx=0, pady=0)

		b = Ribbon.LabelButton(group, image=tkFlair.icons["ok32"],
				text="Ok",
				compound=TOP,
				command=self.ok,
				background=Ribbon._BACKGROUND)
		tkExtra.Balloon.set(b, "Update Material")
		b.pack(side=LEFT, fill=Y, expand=YES)
		self.properties_widgets.append(b)

		b = Ribbon.LabelButton(group, image=tkFlair.icons["x32"],
				text="Cancel",
				compound=TOP,
				command=self.showMaterial,
				background=Ribbon._BACKGROUND)
		tkExtra.Balloon.set(b, "Cancel changes")
		b.pack(side=LEFT, fill=Y, expand=YES)
		self.properties_widgets.append(b)

	#----------------------------------------------------------------------
	# Create Project page
	#----------------------------------------------------------------------
	def createPage(self):
		FlairRibbon.FlairPage.createPage(self)

		self.list_widgets = []
		self.properties_widgets = []

		# ============ Search ============
		frame = Frame(self.frame)
		frame.pack(side=TOP, fill=BOTH, padx=3, pady=2)

		Label(frame, text="Search:").pack(side=LEFT)
		self.searchText = Entry(frame, background="White")
		self.searchText.pack(side=LEFT, expand=YES, fill=X)
		tkExtra.Balloon.set(self.searchText, "Filter materials with search string")
		self.list_widgets.append(self.searchText)

		Button(frame, text="X", padx=1, pady=1,
				command=self.clearSearch).pack(side=RIGHT)

		# ----
		frame = Frame(self.frame)
		frame.pack(side=TOP, fill=BOTH, padx=3, pady=2)

		# -- Group --
		lblframe = LabelFrame(frame, text="Group",
				foreground="DarkBlue", padx=2, pady=2)
		lblframe.pack(side=LEFT, fill=Y)

		self.groupList = tkExtra.ExListbox(lblframe,
					borderwidth=0,
					width=20,
					selectmode=EXTENDED,
					selectborderwidth=0,
					background="White",
					takefocus=True,
					exportselection=FALSE)
		self.groupList.pack(side=LEFT, expand=YES, fill=BOTH)
		self.list_widgets.append(self.groupList)
		sb = Scrollbar(lblframe, orient=VERTICAL, takefocus=False,
				command=self.groupList.yview)
		sb.pack(side=RIGHT, fill=Y)
		self.groupList.config(yscrollcommand=sb.set)

		# -- Material --
		lblframe = LabelFrame(frame, text="Material List",
				foreground="DarkBlue", padx=2, pady=2)
		lblframe.pack(side=RIGHT, expand=YES, fill=BOTH)

		self.materialList = tkExtra.MultiListbox(lblframe,
					(('Material', 30, None),
					 ('Density', 9, None),
					 ('Stoichiometry', 10, None)))
		self.materialList.pack(expand=YES, fill=BOTH)
		self.list_widgets.append(self.materialList)

		# ========== Properties ==========
		lblframe = LabelFrame(self.frame, text="Material Properties",
				foreground="DarkBlue", padx=2, pady=2)
		lblframe.pack(side=BOTTOM, expand=YES, fill=BOTH)

		# ---
		frame = Frame(lblframe)
		frame.pack(side=TOP, fill=X)
		Label(frame, text="Title:").pack(side=LEFT)
		self.title = Entry(frame, background="White",
				disabledforeground="Black")
		tkExtra.Balloon.set(self.title, "Title of the material")
		self.title.pack(side=RIGHT, expand=YES, fill=X)
		self.properties_widgets.append(self.title)

		# === Specific Info ===
		self.tabPages = tkExtra.TabPageSet(lblframe, pageNames=_PAGE_TABS)
		self.tabPages.changePage()
		self.tabPages.pack(side=BOTTOM, expand=YES, fill=BOTH)

		# ========== Notes ==========
		pageframe = self.tabPages[_PAGE_TABS[0]]
		# ---
		row, col = 0, 0
		Label(pageframe, text="Notes:").grid(row=row, column=col, sticky=W)

		# ---
		col += 1
		Label(pageframe, text="Names:").grid(row=row, column=col, sticky=SW)

		row, col = row+1, 0
		rowspan = 6
		self.notes = tkTextEditor.TextEditor(pageframe)
		self.notes.txt.config(width=40, height=5, background="White")
		self.notes.grid(row=row, rowspan=rowspan,
				column=col, #columnspan=colspan,
				padx=3, sticky=NSEW)
		self.properties_widgets.append(self.notes.txt)

		pageframe.grid_columnconfigure(0, weight=1)
		pageframe.grid_rowconfigure(row+rowspan-1, weight=1)

		# --- Names ---
		col += 1
		self.matNameList = tkExtra.ExListbox(pageframe, borderwidth=0,
					selectmode=EXTENDED,
					selectborderwidth=0,
					background="White",
					takefocus=True,
					relief=SUNKEN,
					disabledforeground="Black",
					height=5,
					width=15,
					exportselection=FALSE)
		self.matNameList.grid(row=row, rowspan=rowspan, column=col,
				padx=3, sticky=NSEW)
		self.properties_widgets.append(self.matNameList)

		col += 1
		sb = Scrollbar(pageframe, orient=VERTICAL, takefocus=False,
				command=self.matNameList.yview)
		sb.grid(row=row, rowspan=rowspan, column=col, sticky=NSEW)
		#self.properties_widgets.append(sb)
		self.matNameList.config(yscrollcommand=sb.set)

		# --
		col += 1
		b = Button(pageframe, image=tkFlair.icons["rename"],
				command=self.renameName)
		tkExtra.Balloon.set(b, "Rename")
		b.grid(row=row, column=col, sticky=EW)
		self.properties_widgets.append(b)

		b = Button(pageframe, image=tkFlair.icons["add"],
				command=self.addName)
		tkExtra.Balloon.set(b, "Add name")
		b.grid(row=row+1, column=col, sticky=EW)
		self.properties_widgets.append(b)

		b = Button(pageframe, image=tkFlair.icons["del"],
				command=self.deleteName)
		tkExtra.Balloon.set(b, "Delete name")
		b.grid(row=row+2, column=col, sticky=EW)
		self.properties_widgets.append(b)

		b = Button(pageframe, image=tkFlair.icons["up"],
				command=self.matNameList.moveUp)
		tkExtra.Balloon.set(b, "Move up")
		b.grid(row=row+3, column=col, sticky=EW)
		self.properties_widgets.append(b)

		b = Button(pageframe, image=tkFlair.icons["down"],
				command=self.matNameList.moveDown)
		tkExtra.Balloon.set(b, "Move down")
		b.grid(row=row+4, column=col, sticky=EW)
		self.properties_widgets.append(b)

		# ========== Stoichiometry ==========
		pageframe = self.tabPages[_PAGE_TABS[1]]
		# ---
		row, col = 0, 0
		Label(pageframe, text="Composition:").grid(
				row=row, column=col, sticky=E)
		col += 1
		self.compositionType = tkExtra.Combobox(pageframe, width=7,
				command=self.compositionChanged)
		self.compositionType.fill(_COMPOSITION)
		self.compositionType.config(disabledforeground='Black')
		tkExtra.Balloon.set(self.compositionType, "Type of mixture")
		self.compositionType.grid(row=row, column=col, sticky=W)
		self.properties_widgets.append(self.compositionType)

		# ---
		col += 1
		Label(pageframe, text="State:").grid(
				row=row, column=col, sticky=E)
		col += 1
		self.stateType = tkExtra.Combobox(pageframe, width=7)
		self.stateType.fill(_STATE)
		self.stateType.set(_STATE[0])
		self.stateType.config(disabledforeground='Black')
		tkExtra.Balloon.set(self.stateType, "Material state")
		self.stateType.grid(row=row, column=col, sticky=W)
		self.properties_widgets.append(self.stateType)

		# ---
		col += 1
		Label(pageframe, text="Density:").grid(
				row=row, column=col, sticky=E)
		col += 1
		self.density = tkExtra.FloatEntry(pageframe, background="White",
				disabledforeground="Black")
		tkExtra.Balloon.set(self.density, "Density of material in g/cm^3")
		self.density.grid(row=row, column=col, sticky=W)
		self.properties_widgets.append(self.density)

		# ---
		col += 2
		Label(pageframe, text="Group:").grid(row=row, column=col, sticky=SW)

		# --- Stoichiometry ---
		row, col = row+1, 0
		colspan = 6
		self.stoichiometryList = tkExtra.MultiListbox(pageframe,
					(('Z', 3, None),
					 ('A', 4, None),
					 ('El', 4, None),
					 ('Name', 12, None),
					 ('Frac', 10, None)),
					 height=5)
		self.stoichiometryList.config(disabledforeground='Black')
		self.stoichiometryList.grid(row=row, rowspan=4,
				column=col, columnspan=colspan, sticky=NSEW)
		self.properties_widgets.append(self.stoichiometryList)

		# --
		col += colspan
		b = Button(pageframe, image=tkFlair.icons["edit"],
				command=self.editElement)
		tkExtra.Balloon.set(b, "Edit element")
		b.grid(row=row, column=col, sticky=EW)
		self.properties_widgets.append(b)

		b = Button(pageframe, image=tkFlair.icons["add"],
				command=self.addElement)
		tkExtra.Balloon.set(b, "Add element")
		b.grid(row=row+1, column=col, sticky=EW)
		self.properties_widgets.append(b)

		b = Button(pageframe, image=tkFlair.icons["del"],
				command=self.deleteElement)
		tkExtra.Balloon.set(b, "Delete element")
		b.grid(row=row+2, column=col, sticky=EW)
		self.properties_widgets.append(b)

		# --- Group ---
		col += 1
		rowspan=6
		self.matGroupList = tkExtra.ExListbox(pageframe, borderwidth=0,
					selectmode=MULTIPLE,
					selectborderwidth=0,
					background="White",
					takefocus=True,
					relief=SUNKEN,
					width=15,
					disabledforeground="Black",
					exportselection=FALSE)
		self.matGroupList.grid(row=row, rowspan=rowspan, column=col,
				padx=3, sticky=NSEW)
		self.properties_widgets.append(self.matGroupList)
		col += 1
		sb = Scrollbar(pageframe, orient=VERTICAL, takefocus=False,
				command=self.matGroupList.yview)
		sb.grid(row=row, rowspan=rowspan, column=col, sticky=NSEW)
		#self.properties_widgets.append(sb)
		self.matGroupList.config(yscrollcommand=sb.set)

		# ---
		pageframe.grid_rowconfigure(row+3, weight=1)
		pageframe.grid_columnconfigure(2, weight=1)
		pageframe.grid_columnconfigure(4, weight=1)

		# --- Frame grid weights ---
#		mainframe.grid_columnconfigure(0, weight=1)
#		mainframe.grid_columnconfigure(1, weight=5)
#		mainframe.grid_rowconfigure(0, weight=1)
#		mainframe.grid_rowconfigure(1, weight=3)

		# --- Bindings ---
		self.searchText.bind('<Return>',         self.fillMaterials)
		self.searchText.bind('<KP_Enter>',       self.fillMaterials)
		self.groupList.bind('<<ListboxSelect>>', self.fillMaterials)
		#self.materialList.bindList('<ButtonRelease-1>', self.showMaterial)
		self.materialList.bind('<<ListboxSelect>>', self.showMaterial)

		# Doesn't work!
		self.materialList.bindList('<Double-1>', self.insert2Input)
		self.materialList.bindList('<Return>',   self.insert2Input)
		self.materialList.bindList('<KP_Enter>', self.insert2Input)
		self.matNameList.bind('<Double-1>',      self.renameName)
		self.matNameList.bind('<Control-Key-E>', self.renameName)

		# --- Menus ---
		self.materialList.setPopupMenu([
				("Insert to Input", 0, self.insert2Input,   tkFlair.icons["input"]),
				None,
				("Edit",   0, self.edit,   tkFlair.icons["edit"]),
				("Add",    0, self.add,    tkFlair.icons["add"]),
				("Delete", 0, self.delete, tkFlair.icons["x"])])

		self.matNameList.setPopupMenu(
				[("Rename", 0, self.renameName,  tkFlair.icons["rename"]),
				 ("Add",    0, self.addName,     tkFlair.icons["add"]),
				 ("Delete", 0, self.deleteName,  tkFlair.icons["x"])])

		self.stoichiometryList.setPopupMenu(
				[("Edit",   0, self.editElement,   tkFlair.icons["edit"]),
				 ("Add",    0, self.addElement,    tkFlair.icons["add"]),
				 ("Delete", 0, self.deleteElement, tkFlair.icons["x"])])

		self.stoichiometryList.bindList('<Double-1>', self.editElement)
		self.stoichiometryList.bindLeftRight()
		self.stoichiometryList.bindList('<Return>',   self.editElement)
		self.stoichiometryList.bindList('<KP_Enter>', self.editElement)

		# ========== Properties ==========
		pageframe = self.tabPages[_PAGE_TABS[2]]

		# ---
		row, col = 0, 0
		#Label(pageframe, text="State:").grid(
		#		row=row, column=col, sticky=E)
		#col += 1
		#self.stateType = StringVar()
		#self.stateType.set(_STATE[0])
		#o = tkExtra.ExOptionMenu(pageframe, self.stateType,
		#			*_STATE)
		#o.config(disabledforeground='Black')
		#tkExtra.Balloon.set(o, "Material state")
		#o.grid(row=row, column=col, sticky=W)
		#self.properties_widgets.append(o)
#
#		# ---
#		row, col = row+1, 0
		Label(pageframe, text="Pressure (atm):").grid(
				row=row, column=col, sticky=E)
		col += 1
		self.pressure = tkExtra.FloatEntry(pageframe, background="White",
				disabledforeground="Black")
		tkExtra.Balloon.set(self.pressure, "Gas pressure in atmospheres")
		self.pressure.grid(row=row, column=col, sticky=W)
		self.properties_widgets.append(self.pressure)

		# ---
		row, col = row+1, 0
		Label(pageframe, text="RHOR Factor:").grid(
				row=row, column=col, sticky=E)
		col += 1
		self.rhor = tkExtra.FloatEntry(pageframe, background="White",
				disabledforeground="Black")
		tkExtra.Balloon.set(self.rhor,
			"Factor to multiply the density of the material, "
			"when calculating density effect parameters")
		self.rhor.grid(row=row, column=col, sticky=W)
		self.properties_widgets.append(self.rhor)

		# ---
		row, col = row+1, 0
		Label(pageframe, text="Ionization Potential (eV):").grid(
				row=row, column=col, sticky=E)
		col += 1
		self.ionization = tkExtra.FloatEntry(pageframe, background="White",
				disabledforeground="Black")
		tkExtra.Balloon.set(self.ionization,
			"Average Ionization potential to be used for dE/dx")
		self.ionization.grid(row=row, column=col, sticky=W)
		self.properties_widgets.append(self.ionization)

		# ---
		row, col = row+1, 0
		Label(pageframe, text="DPA Atomic Displacement Energy (eV):").grid(
				row=row, column=col, sticky=E)
		col += 1
		self.dpa_threshold = tkExtra.FloatEntry(pageframe, background="White",
				disabledforeground="Black")
		tkExtra.Balloon.set(self.dpa_threshold,
			"DPA Displacement energy threshold averaged " \
			"over all crystallographic directions")
		self.dpa_threshold.grid(row=row, column=col, sticky=W)
		self.properties_widgets.append(self.dpa_threshold)

		# ---
		row, col = row+1, 0
		Label(pageframe, text="Density scaling factor for dE/dx:").grid(
				row=row, column=col, sticky=E)
		col += 1
		self.xrho_dEdx = tkExtra.FloatEntry(pageframe, background="White",
				disabledforeground="Black")
		tkExtra.Balloon.set(self.xrho_dEdx,
			"density scaling factor for charged particle ionisation " \
			"processes (dE/dx, delta ray production, Moller and Bhabha " \
			"scattering")
		self.xrho_dEdx.grid(row=row, column=col, sticky=W)
		self.properties_widgets.append(self.xrho_dEdx)

		# ---
		row, col = row+1, 0
		Label(pageframe, text="Density scaling factor for other:").grid(
				row=row, column=col, sticky=E)
		col += 1
		self.xrho = tkExtra.FloatEntry(pageframe, background="White",
				disabledforeground="Black")
		tkExtra.Balloon.set(self.xrho,
			"density scaling factor for all other processes")
		self.xrho.grid(row=row, column=col, sticky=W)
		self.properties_widgets.append(self.xrho)

		# ---
		row = 0
		col += 1
		rowspan = 6
		frame = LabelFrame(pageframe, text="Sternheimer parameters",
				foreground="DarkBlue", padx=2, pady=2)
		frame.grid(row=row, rowspan=rowspan, column=col, sticky=NSEW)

		# ---
		row, col = 0,0
		Label(frame, text="Cbar:").grid(
				row=row, column=col, sticky=E)
		col += 1
		self.sternheimer_cbar = tkExtra.FloatEntry(frame, background="White",
				disabledforeground="Black")
		tkExtra.Balloon.set(self.sternheimer_cbar, "Sternheimer Cbar parameter")
		self.sternheimer_cbar.grid(row=row, column=col, sticky=W)
		self.properties_widgets.append(self.sternheimer_cbar)

		# ---
		row, col = row+1, 0
		Label(frame, text="X0:").grid(
				row=row, column=col, sticky=E)
		col += 1
		self.sternheimer_x0 = tkExtra.FloatEntry(frame, background="White",
				disabledforeground="Black")
		tkExtra.Balloon.set(self.sternheimer_x0, "Sternheimer X0 parameter")
		self.sternheimer_x0.grid(row=row, column=col, sticky=W)
		self.properties_widgets.append(self.sternheimer_x0)

		# ---
		row, col = row+1, 0
		Label(frame, text="X1:").grid(
				row=row, column=col, sticky=E)
		col += 1
		self.sternheimer_x1 = tkExtra.FloatEntry(frame, background="White",
				disabledforeground="Black")
		tkExtra.Balloon.set(self.sternheimer_x1, "Sternheimer X1 parameter")
		self.sternheimer_x1.grid(row=row, column=col, sticky=W)
		self.properties_widgets.append(self.sternheimer_x1)

		# ---
		row, col = row+1, 0
		Label(frame, text="a:").grid(
				row=row, column=col, sticky=E)
		col += 1
		self.sternheimer_a = tkExtra.FloatEntry(frame, background="White",
				disabledforeground="Black")
		tkExtra.Balloon.set(self.sternheimer_a, "Sternheimer a parameter")
		self.sternheimer_a.grid(row=row, column=col, sticky=W)
		self.properties_widgets.append(self.sternheimer_a)

		# ---
		row, col = row+1, 0
		Label(frame, text="m:").grid(
				row=row, column=col, sticky=E)
		col += 1
		self.sternheimer_m = tkExtra.FloatEntry(frame, background="White",
				disabledforeground="Black")
		tkExtra.Balloon.set(self.sternheimer_m, "Sternheimer m (SK) parameter")
		self.sternheimer_m.grid(row=row, column=col, sticky=W)
		self.properties_widgets.append(self.sternheimer_m)

		# ---
		row, col = row+1, 0
		Label(frame, text="\u03b40:").grid(
				row=row, column=col, sticky=E)
		col += 1
		self.sternheimer_d0 = tkExtra.FloatEntry(frame, background="White",
				disabledforeground="Black")
		tkExtra.Balloon.set(self.sternheimer_d0, "Sternheimer Delta0 parameter")
		self.sternheimer_d0.grid(row=row, column=col, sticky=W)
		self.properties_widgets.append(self.sternheimer_d0)

		# ================================
		tkExtra.multiConfig(self.properties_widgets, state=DISABLED)
		self.editing = False
		self.shown   = -1
		self.compositionValue = ""
		self.refresh()

	# ----------------------------------------------------------------------
	# Show Frame
	# ----------------------------------------------------------------------
	def showFrame(self):
		tkFlair.TitleFrame.showFrame(self)
		self.refresh()

	# ----------------------------------------------------------------------
	def refresh(self):
		FlairRibbon.FlairPage.refresh(self)
		self.fillGroups()
		self.groupList.selection_set(0, END)
		self.fillMaterials()

	# ----------------------------------------------------------------------
	def fillGroups(self, event=None):
		# Fill groups
		self.groupList.delete(0, END)

		glist = groups[:]
		glist.sort()
		for g in glist:
			self.groupList.insert(END, g)

	# ----------------------------------------------------------------------
	def clearSearch(self):
		self.searchText.delete(0,END)
		self.fillMaterials()

	# ----------------------------------------------------------------------
	def fillMaterials(self, event=None):
		# Selected groups
		selgroup = []
		for sel in map(int, self.groupList.curselection()):
			group = self.groupList.get(sel)
			gid = groups.index(group) + 1
			selgroup.append(gid)

		# Search string
		search_string = self.searchText.get().strip().upper().split()
		if len(search_string) == 0: search_string = None

		# Fill list
		self.materialList.delete(0, END)
		savesort = self.materialList.saveSort()

		for mat in list(materials.values()):
			# Find intersection of groups
			for mg in mat.group:
				if mg in selgroup: break
			else:
				continue

			# Check also search pattern
			if search_string:
				if Utils.search(mat.title.upper(), search_string):
					pass
				elif Utils.search(mat.stoichiometry.upper(), search_string):
					pass
				elif Utils.search(mat.notes.upper(), search_string):
					pass
				elif Utils.search(mat.name.upper(), search_string):
					pass
				else:
					continue

			# Add material
			self.materialList.insert(END,
					(mat.title, mat.density, mat.stoichiometry))

		self.materialList.restoreSort(savesort)
		self.materialList.selection_clear(0, END)
		self.materialList.selection_set(0)
		self.showMaterial()

	# ----------------------------------------------------------------------
	def get(self, event=None):
		if self.page is None: return
		if self.editing and self.materialChanged():
			# -----------------------------------
			# Ask user to update the old material
			# -----------------------------------
			if tkFlair.askyesno("Update material: "+self.material.title,
			   "Material \"%s\" has changed.\n"
			   "Do you want to update the material?"%(self.material.title),
			   parent=self.page):
				self.updateMaterial()

		self.editing = False
		tkExtra.multiConfig(self.list_widgets, state=NORMAL)
		tkExtra.multiConfig(self.properties_widgets, state=DISABLED)

	# ----------------------------------------------------------------------
	def showMaterial(self, event=None):
		self.get()

		# Empty material info
		tkExtra.multiConfig(self.properties_widgets, state=NORMAL)
		self.title.delete(0, END)
		self.density.delete(0, END)
		self.notes.txt.delete(1.0, END)
		self.stoichiometryList.delete(0, END)
		self.matGroupList.delete(0, END)
		self.matNameList.delete(0, END)
		self.pressure.delete(0, END)
		self.rhor.delete(0, END)
		self.ionization.delete(0, END)
		self.dpa_threshold.delete(0, END)
		self.xrho_dEdx.delete(0, END)
		self.xrho.delete(0, END)
		self.sternheimer_cbar.delete(0, END)
		self.sternheimer_x0.delete(0, END)
		self.sternheimer_x1.delete(0, END)
		self.sternheimer_a.delete(0, END)
		self.sternheimer_m.delete(0, END)
		self.sternheimer_d0.delete(0, END)
		self.compositionType.set("")
		self.stateType.set("")

		# Show selected one
		self.shown = -1
		sel = self.materialList.curselection()
		if len(sel) != 1:
			tkExtra.multiConfig(self.properties_widgets, state=DISABLED)
			return

		# Search for material in database
		self.shown = int(sel[0])
		(title, density, stoichiometry) = self.materialList.get(self.shown)
		mat = materials[title]
		self.material = mat

		# Fill values
		self.title.insert(0, mat.title)
		self.density.insert(0, mat.density)
		self.notes.txt.insert(0.0, mat.notes)
		self.compositionType.set(_COMPOSITION[mat.composition])
		self.compositionValue = ""	# Empty
		self.stateType.set(_STATE[mat.state])

		self.pressure.insert(         0, mat.pressure)
		self.rhor.insert(             0, mat.rhor)
		self.ionization.insert(       0, mat.ionization)
		self.dpa_threshold.insert(    0, mat.dpa_threshold)
		self.xrho_dEdx.insert(        0, mat.xrho_dEdx)
		self.xrho.insert(             0, mat.xrho)

		self.sternheimer_cbar.insert( 0, mat.s_cbar)
		self.sternheimer_x0.insert(   0, mat.s_x0)
		self.sternheimer_x1.insert(   0, mat.s_x1)
		self.sternheimer_a.insert(    0, mat.s_a)
		self.sternheimer_m.insert(    0, mat.s_m)
		self.sternheimer_d0.insert(   0, mat.s_d0)

		# Fill Stoichiometry
		for Z, A, f in mat.elements:
			self.insertElement(END, Z, A, f)

		# Fill groups
		for g in mat.group:
			self.matGroupList.insert(END, groups[g-1])

		# Fill names
		for name in mat.name.split():
			self.matNameList.insert(END, name)

		tkExtra.multiConfig(self.properties_widgets, state=DISABLED)

	# ----------------------------------------------------------------------
	def edit(self, event=None):
		if self.shown<0: return

		self.editing = True

		tkExtra.multiConfig(self.list_widgets, state=DISABLED)
		tkExtra.multiConfig(self.properties_widgets, state=NORMAL)

		# Remember composition type
		self.compositionValue = _COMPOSITION[self.material.composition]

		# Change group list to BROWSE
		self.title.focus_set()
		self.matGroupList.delete(0, END)
		glist = groups[:]
		glist.sort()
		row = 0
		for g in glist:
			self.matGroupList.insert(END, g)
			gi = groups.index(g)+1
			if gi in self.material.group:
				self.matGroupList.selection_set(row)
			row += 1

	# ----------------------------------------------------------------------
	def compositionChanged(self):
		if self.compositionValue == "": return
		new = self.compositionType.get()
		if new == "" or new == self.compositionValue: return

		old = self.compositionValue
		self.compositionValue = new

		# XXX Ask the user before?
		if not tkFlair.askyesno("Composition type changed",
			"Do you want to update the elements fractions\n"
			"from \"%s\" to \"%s\"?"%(old, new), parent=self.page):
			return

		# Convert everything to mass
		flist = []
		for i in range(self.stoichiometryList.size()):
			Z, A, S, N, f = self.stoichiometryList.get(i)
			element = PeriodicTable.element(Z)
			if A == "": A = element.Amass()
			mf = 0.0
			if   old == "atom":
				try:
					mf = f*float(A)
				except:
					pass
			elif old == "volume":
				try:
					mf = f*element.density()
				except:
					pass
			else:
				mf = f

			# And then to the new type
			if   new == "atom":
				try:
					f = mf/float(A)
				except:
					pass
			elif new == "volume":
				try:
					f = mf/element.density()
				except:
					pass
			else:
				f = mf
			flist.append(f)

		# Insert items
		sel = list(map(int, self.stoichiometryList.curselection()))
		listbox = self.stoichiometryList.lists[-1]
		listbox.delete(0, END)
		i = 0
		for f in flist:
			listbox.insert(END, f)
			if i in sel: listbox.selection_set(i)
			i += 1

	# ----------------------------------------------------------------------
	def ok(self):
		self.updateMaterial()
		self.showMaterial()

	# ----------------------------------------------------------------------
	# Check if material description has changed
	# ----------------------------------------------------------------------
	def materialChanged(self, event=None):
		mat = self.material

		if mat.title != self.title.get().strip(): return True

		if mat.notes != self.notes.txt.get(1.0, END).strip(): return True

		name = " ".join(list(self.matNameList.get(0,END)))
		if mat.name  != name: return True

		composition = _COMPOSITION.index(self.compositionType.get())
		if mat.composition != composition: return True

		# find groups
		group = []
		for g in self.matGroupList.curselection():
			gname = self.matGroupList.get(g)
			group.append(groups.index(gname)+1)
		group.sort()

		if " ".join(map(str,mat.group)) != " ".join(map(str,group)):
			return True

		# Do the checksum
		chksum = 0.0
		try: chksum += float(self.density.get())
		except: pass

		try: chksum += float(self.pressure.get())
		except: pass
		try: chksum += float(self.rhor.get())
		except: pass
		try: chksum += float(self.ionization.get())
		except: pass
		try: chksum += float(self.dpa_threshold.get())
		except: pass
		try: chksum += float(self.xrho_dEdx.get())
		except: pass
		try: chksum += float(self.xrho.get())
		except: pass

		try: chksum += float(self.sternheimer_cbar.get())
		except: pass
		try: chksum += float(self.sternheimer_x0.get())
		except: pass
		try: chksum += float(self.sternheimer_x1.get())
		except: pass
		try: chksum += float(self.sternheimer_a.get())
		except: pass
		try: chksum += float(self.sternheimer_m.get())
		except: pass
		try: chksum += float(self.sternheimer_d0.get())
		except: pass

		for row in range(self.stoichiometryList.size()):
			Z,A,symbol,name,f = self.stoichiometryList.get(row)
			if A=="": A=0
			try: chksum += float(Z)
			except: pass
			try: chksum += float(A)
			except: pass
			try: chksum += float(f)
			except: pass

		if abs(chksum - mat.chksum) > 1e-10: return True

		return False

	# ----------------------------------------------------------------------
	def updateMaterial(self):
		global materials, changed
		mat = self.material

		newtitle = self.title.get()
		composition = _COMPOSITION.index(self.compositionType.get())
		if newtitle != mat.title:
			del materials[mat.title]
			materials[newtitle] = mat
			mat.title = newtitle

		mat.name  = " ".join(list(self.matNameList.get(0,END)))
		mat.notes = self.notes.txt.get(1.0, END).strip()

		try: mat.density     = float(self.density.get())
		except: mat.density  = 0.0

		try: mat.pressure    = float(self.pressure.get())
		except: mat.pressure = ""
		try: mat.rhor        = float(self.rhor.get())
		except: mat.rhor     = ""
		try: mat.ionization  = float(self.ionization.get())
		except: mat.ionization = ""
		try: mat.dpa_threshold = float(self.dpa_threshold.get())
		except: mat.dpa_threshold = ""
		try: mat.xrho_dEdx = float(self.xrho_dEdx.get())
		except: mat.xrho_dEdx = ""
		try: mat.xrho= float(self.xrho.get())
		except: mat.xrho= ""

		try: mat.s_cbar  = float(self.sternheimer_cbar.get())
		except: mat.s_cbar = ""
		try: mat.s_x0    = float(self.sternheimer_x0.get())
		except: mat.s_x0 = ""
		try: mat.s_x1    = float(self.sternheimer_x1.get())
		except: mat.s_x1 = ""
		try: mat.s_a     = float(self.sternheimer_a.get())
		except: mat.s_a  = ""
		try: mat.s_m     = float(self.sternheimer_m.get())
		except: mat.s_m  = ""
		try: mat.s_d0    = float(self.sternheimer_d0.get())
		except: mat.s_d0 = ""

		mat.group = []
		for g in self.matGroupList.curselection():
			gname = self.matGroupList.get(g)
			mat.group.append(groups.index(gname)+1)
		mat.group.sort()
		mat.composition = _COMPOSITION.index(self.compositionType.get())
		mat.state       = _STATE.index(self.stateType.get())

		mat.elements = []
		for row in range(self.stoichiometryList.size()):
			Z,A,symbol,name,f = self.stoichiometryList.get(row)
			if A=="": A=0
			Z = int(Z)
			A = int(A)
			f = float(f)
			mat.elements.append((Z,A,f))

		mat.make()
		mat.setModified()
		changed = True

		# Update listbox
		self.materialList.config(state=NORMAL)
		active = self.materialList.index(ACTIVE)
		self.materialList.delete(self.shown)
		self.materialList.insert(self.shown,
				(mat.title, mat.density, mat.stoichiometry))
		if active == self.shown:
			self.materialList.selection_set(self.shown)
			self.materialList.activate(self.shown)

	# ----------------------------------------------------------------------
	def add(self, event=None):
		global materials, changed
		if self.editing: return
		changed = True
		mat = Material("New Material")
		mat.name = "NewMat"
		mat.group = [1]		# User
		mat.make()
		materials[mat.title] = mat
		self.materialList.insert(END,
			(mat.title, mat.density, mat.stoichiometry))
		self.materialList.selection_clear(0, END)
		self.materialList.selection_set(END)
		self.materialList.activate(END)
		self.materialList.see(END)
		self.showMaterial()
		self.edit()

	# ----------------------------------------------------------------------
	def delete(self, event=None):
		if self.editing: return
		global materials, changed
		if not tkFlair.askyesno("Delete materials",
				"Do you want to delete the selected materials?", parent=self.page):
			return
		changed = True
		sel = list(self.materialList.curselection())
		sel.reverse()
		for i in sel:
			mat = self.materialList.get(i)
			del materials[mat[0]]
			self.materialList.delete(i)

	# ----------------------------------------------------------------------
	def editElement(self, event=None):
		if not self.editing: return
		listbox = self.page.focus_get()
		if listbox is None: return

		try: active = listbox.index(ACTIVE)
		except: return
		values = self.stoichiometryList.get(ACTIVE)
		lid    = self.stoichiometryList.lists.index(listbox)

		# Edit and change value
		if   lid < 2:
			# Edit Z or A
			edit = tkExtra.InPlaceInteger(listbox)
			if edit.value is None or lid==1:
				return
			try:
				Z = int(edit.value)
			except:
				return
			element = PeriodicTable.element(Z)

		elif lid < 4:
			# Edit Symbol or Name
			edit = InPlaceElement(listbox)
			if edit.value is None:
				return
			element = PeriodicTable.elementBySymbol(edit.value)
			if element is None:
				element = PeriodicTable.elementByName(edit.value)
			if element:
				Z = element.Z
			else:
				return

		else:
			# Edit fraction
			edit = tkExtra.InPlaceFloat(listbox)
			return

		# Update active row
		self.stoichiometryList.set(active,
		(Z, values[1], element.symbol, element.name, values[4]))
		self.stoichiometryList.selection_set(active)
		self.stoichiometryList.activate(active)

	# ----------------------------------------------------------------------
	def addElement(self, event=None):
		self.insertElement(END, 1, 0, 0.0)
		self.stoichiometryList.selection_clear(0,END)
		self.stoichiometryList.selection_set(END)
		self.stoichiometryList.activate(END)
		self.stoichiometryList.focus_set()
		self.editElement()

	# ----------------------------------------------------------------------
	def insertElement(self, pos, Z, A, f):
		if A==0: A=""
		element = PeriodicTable.element(Z)
		self.stoichiometryList.insert(pos,
			(Z, A, element.symbol, element.name, f))

	# ----------------------------------------------------------------------
	def deleteElement(self, event=None):
		sel = list(self.stoichiometryList.curselection())
		if len(sel) == 0: return
		sel.reverse()
		for i in sel:
			self.stoichiometryList.delete(i)

	# ----------------------------------------------------------------------
	def renameName(self, event=None):
		if self.matNameList["state"] == NORMAL:
			edit = tkExtra.InPlaceEdit(self.matNameList)

	# ----------------------------------------------------------------------
	def addName(self, event=None):
		self.matNameList.insert(END,"")
		self.matNameList.selection_clear(0,END)
		self.matNameList.selection_set(END)
		self.matNameList.activate(END)
		self.renameName()
		if self.matNameList.get(ACTIVE).strip() == "":
			self.matNameList.delete(ACTIVE)
			self.matNameList.selection_set(END)
			self.matNameList.activate(END)

	# ----------------------------------------------------------------------
	def deleteName(self, event=None):
		sel = list(self.matNameList.curselection())
		if len(sel) == 0: return
		sel.reverse()
		for i in sel:
			self.matNameList.delete(i)

	# ----------------------------------------------------------------------
	# Insert selected materials to input
	# ----------------------------------------------------------------------
	def insert2Input(self, event=None):
		# Insert selected materials
		mats = []
		for sel in self.materialList.curselection():
			mats.append(self.materialList.get(sel)[0])
		matnames, elements = insert2Input(self.flair, mats)
		if matnames: self.setModified()

	# ----------------------------------------------------------------------
	# Import materials from input
	# ----------------------------------------------------------------------
	def importFromInput(self, event=None):
		global materials, changed

		materialDict = self.project.input.materialDict()
		material = self.project.input["MATERIAL"]
		compound = self.project.input["COMPOUND"]

		errorstr = ""
		imported = []

		for name in SelectInputMaterialsDialog(self.page, self.project.input).show():
			# Find material card
			for card in material:
				if card.sdum() == name: break
			else:
				errorstr += "Material %r card not found\n"%(name)
				continue

			mat = Material(name)
			mat.name    = name
			mat.group   = [1]		# User
			mat.notes   = card.comment()
			mat.density = card.numWhat(3)

			for card in compound:
				if card.sdum() != name: continue
				for i in range(1,6,2):
					f = card.numWhat(i)
					iso  = card.what(i+1)
					if not iso: continue
					if f>0.0:
						self.composition = 0
					elif iso[0]=="-":
						iso = iso[1:]
						self.composition = 2
					else:
						f = -f
						self.composition = 1

					# Find isotope
					ciso = materialDict.get(iso)
					if not ciso:
						errorstr += "Material %r Isotope %r not found\n"%(name,iso)
						continue
					Z = ciso.intWhat(1)
					if Z==0:
						errorstr += "Material %r cannot handle compound with compounds %r\n"%(name,iso)
					A = ciso.intWhat(6)
					mat.elements.append((Z,A,f))
			imported.append(name)

			mat.setModified()
			mat.make()
			materials[mat.title] = mat
			changed = True

		if imported:
			if errorstr:
				self.flair.notify("Imported Materials",
					"%d materials imported with errors\n%s"%(len(imported),errorstr),
					tkFlair.NOTIFY_WARNING)
			else:
				self.flair.notify("Imported Materials",
					"%d materials imported"%(len(imported)))
		elif errorstr:
			self.flair.notify("Import Errors",
					errorstr,
					tkFlair.NOTIFY_ERROR)

#		self.materialList.insert(END,
#			(mat.title, mat.density, mat.stoichiometry))
#		self.materialList.selection_clear(0, END)
#		self.materialList.selection_set(END)
#		self.materialList.activate(END)
#		self.materialList.see(END)

	# ----------------------------------------------------------------------
	def copy(self, event=None):
		global _clipboard

		# XXX copy material, stoichiometryList etc...
		fw = self.page.focus_get()
		if fw in self.stoichiometryList.lists:
			del _clipboard[:]
			sel = list(map(int, self.stoichiometryList.curselection()))
			for i in sel:
				(Z,A,S,N,f) = self.stoichiometryList.get(i)
				_clipboard.append((Z,A,f))
		else:
			fw.event_generate("<<Copy>>")

	# ----------------------------------------------------------------------
	def cut(self, event=None):
		# XXX cut material, stoichiometryList etc...
		fw = self.page.focus_get()
		if fw in self.stoichiometryList.lists:
			self.copy()
			self.deleteElement()
		else:
			fw.event_generate("<<Cut>>")

	# ----------------------------------------------------------------------
	def paste(self, event=None):
		global _clipboard

		# XXX paste material, stoichiometryList etc...
		fw = self.page.focus_get()
		if fw in self.stoichiometryList.lists:
			pos = END
			sel = self.stoichiometryList.curselection()
			if len(sel)>0: pos = int(sel[-1])+1
			for Z, A, f in _clipboard:
				self.insertElement(pos, Z, A, f)
				if pos != END: pos += 1
			self.stoichiometryList.selection_clear(0,END)
			self.stoichiometryList.set(pos)
			self.stoichiometryList.activate(pos)
			self.stoichiometryList.see(pos)
		else:
			fw.event_generate("<<Paste>>")
			fw.focus_set()

	# ----------------------------------------------------------------------
	def find(self):
		self.notes.findDialog()

	# ----------------------------------------------------------------------
	def replace(self):
		self.notes.replaceDialog()

	# ----------------------------------------------------------------------
	def findNext(self):
		self.notes.findNext()

#===============================================================================
class SelectPropertiesDlg(Toplevel):
	def __init__(self, master, material, prop):
		Toplevel.__init__(self, master)
		self.title("%s properties"%(material))

		l = Label(self, text= \
			"Select the additional properties of "
			"material \"%s\" to be added in the input"%(material),
			foreground="DarkBlue", background="White",
			justify=CENTER)
		l.pack(side=TOP, fill=X)

		self.prop = prop
		self.vars = {}

		v = self._addButton(self,
			"MAT-PROP", "pressure, ionization potential")
		v.pack(side=TOP, fill=X)

		v = self._addButton(self,
			"DPA", "DPA damage threshold")
		v.pack(side=TOP, fill=X)

		v = self._addButton(self,
			"STERNHEI", "Sternheimer parameters")
		v.pack(side=TOP, fill=X)

		v = self._addButton(self,
			"CORRFACT", "Density correction factor")
		v.pack(side=TOP, fill=X)

		frame = Frame(self)
		frame.pack(side=BOTTOM, expand=YES, fill=X)
		self.bind("<Escape>", self.close)
		b = Button(frame, text="Only MATERIAL card", command=self.close)
		b.pack(side=RIGHT)
		b = Button(frame, text="All the above", command=self.insertAll)
		b.pack(side=RIGHT)
		b = Button(frame, text="Selected Only", command=self.insert)
		b.pack(side=RIGHT)

		self.add = []
		x = master.winfo_rootx() + master.winfo_width()/2 -100
		y = master.winfo_rooty() + master.winfo_height()/2 -50
		self.geometry("+%d+%d" %(x,y))
		self.focus_set()
		self.wait_visibility()
		self.grab_set()
		self.wait_window()

	# ----------------------------------------------------------------------
	def _addButton(self, frame, var, text):
		v = IntVar()
		v.set(0)
		c = Checkbutton(frame, text="%s: %s"%(var,text),
			variable=v, anchor=W, justify=LEFT,
			background="White")
		if var not in self.prop:
			c["state"] = DISABLED
		self.vars[var] = v
		return c

	# ----------------------------------------------------------------------
	def insert(self):
		for k,v in list(self.vars.items()):
			if v.get():
				self.add.append(k)
		if not self.add:
			messagebox.showinfo("None selected",
				"No property is selected for insertion in the input",
				parent=self)
			return
		self.destroy()

	# ----------------------------------------------------------------------
	def insertAll(self):
		self.add = self.prop
		self.destroy()

	# ----------------------------------------------------------------------
	def close(self, event=None):
		self.destroy()

#===============================================================================
class SelectInputMaterialsDialog(Toplevel):
	def __init__(self, master, input):
		Toplevel.__init__(self, master)
		self.transient(master)

		# Listbox
		frame = Frame(self)
		frame.pack(expand=YES,fill=BOTH)

		self.listbox = tkExtra.ExListbox(frame, borderwidth=0,
			selectmode=MULTIPLE,
			selectborderwidth=0,
			background="White",
			takefocus=True,
			exportselection=FALSE)
		self.listbox.pack(side=LEFT, expand=YES, fill=BOTH)

		sb = Scrollbar(frame, orient=VERTICAL, command=self.listbox.yview)
		sb.pack(side=RIGHT, fill=Y)
		self.listbox.config(yscrollcommand=sb.set)

		frame = Frame(self)
		frame.pack(side=BOTTOM)

		b = Button(frame,text="Cancel",
				image=tkFlair.icons["x"],
				compound=LEFT,
				takefocus=FALSE,
				command=self.cancel)
		b.pack(side=RIGHT)
		b = Button(frame,text="All",
				image=tkFlair.icons["all"],
				compound=LEFT,
				takefocus=FALSE,
				command=self.all)
		b.pack(side=RIGHT)
		b = Button(frame,text="Ok",
				image=tkFlair.icons["ok"],
				compound=LEFT,
				takefocus=FALSE,
				command=self.ok)
		b.pack(side=RIGHT)
		self.bind("<Escape>", self.cancel)
		self.bind("<Control-Key-a>", self.all)
		self.protocol("WM_DELETE_WINDOW", self.cancel)
		self.title("Import Materials")

		# Fill template box
		mats = {}
		for card in input["COMPOUND"]:
			mats[card.sdum()] = True
		mats = list(mats.keys())
		mats.sort()
		for m in mats:
			self.listbox.insert(END, m)
		self.selection = []

	# ----------------------------------------------------------------------
	def show(self):
		self.focus_set()
		self.listbox.focus_set()
		self.wait_visibility()
		self.grab_set()
		self.wait_window()
		return self.selection

	# ----------------------------------------------------------------------
	def all(self, event=None):
		self.listbox.selection_set(0,END)
		self.ok()

	# ----------------------------------------------------------------------
	def ok(self, event=None):
		self.update_idletasks()
		self.selection = [self.listbox.get(x) for x in self.listbox.curselection()]
		self.destroy()

	# ----------------------------------------------------------------------
	def cancel(self, event=None):
		self.destroy()

#===============================================================================
class InPlaceElement(tkExtra.InPlaceEdit):
	# ----------------------------------------------------------------------
	def createWidget(self):
		self.edit = Entry(self.frame, width=10, background="White")
		self.edit.pack(side=LEFT, expand=YES, fill=BOTH)
		b = Button(self.frame, image=tkFlair.icons["toi"],
			padx=0, pady=0, command=self.showTable)
		b.pack(side=RIGHT)
		self.edit.focus_set()

	# ----------------------------------------------------------------------
	def showTable(self):
		self.frame.unbind("<FocusOut>")
		element = PeriodicTable.show(None, False, True)
		self.frame.bind("<FocusOut>", self.cancel)
		self.edit.focus_set()
		if element is not None:
			self.edit.delete(0,END)
			self.edit.insert(0, element.symbol)

#-------------------------------------------------------------------------------
# Initialize
#-------------------------------------------------------------------------------
def init(filenames):
	global materials, groups

	# read card database and prepare CardInfo classes
	matini = configparser.RawConfigParser()
	if matini.read(filenames) is []:
		raise Exception("Cannot open file "+filenames)

	# Read groups...
	i = 1
	try:
		while True:
			g = matini.get("groups","group.%d"%(i))
			groups.append(g)
			i += 1
	except configparser.NoOptionError:
		pass

	# Read Materials
	for material in matini.sections():
		if material == "groups": continue
		mat = Material(material)

		mat.name  = matini.get(material, "name")
#		except: pass

		try:
			mat.group = list(map(int, matini.get(material, "groups").split()))
			mat.group.sort()
		except: pass

		try: mat.density = matini.getfloat(material, "density")
		except: pass

		try: mat.composition = _COMPOSITION.index(matini.get(material, "composition"))
		except: pass

		try:
			notes = matini.get(material, "notes").decode("utf-8")
			mat.notes = string.replace(notes[1:],"\n$","\n").strip()
		except: pass

		try: mat.modified = matini.get(material, "modified")
		except: pass

		try: mat.state  = _STATE.index(matini.get(material, "state"))
		except: pass
		try: mat.pressure = matini.getfloat(material, "pressure")
		except: pass
		try: mat.rhor     = matini.getfloat(material, "rhor")
		except: pass
		try: mat.ionization = matini.getfloat(material, "ionization")
		except: pass
		try: mat.dpa_threshold = matini.getfloat(material, "dpa")
		except: pass
		try: mat.xrho_dEdx = matini.getfloat(material, "xrho_dEdx")
		except: pass
		try: mat.xrho      = matini.getfloat(material, "xrho")
		except: pass

		try: mat.s_cbar = matini.getfloat(material, "s_cbar")
		except: pass
		try: mat.s_x0   = matini.getfloat(material, "s_x0")
		except: pass
		try: mat.s_x1   = matini.getfloat(material, "s_x1")
		except: pass
		try: mat.s_a    = matini.getfloat(material, "s_a")
		except: pass
		try: mat.s_m    = matini.getfloat(material, "s_m")
		except: pass
		try: mat.s_d0   = matini.getfloat(material, "s_d0")
		except: pass

		n = matini.getint(material, "elements")
		for i in range(1, n+1):
			Z = matini.getint(material, "z.%d"%(i))
			try: A = matini.getint(material, "a.%d"%(i))
			except: A = 0
			f = matini.getfloat(material, "f.%d"%(i))
			mat.elements.append((Z, A, f))
		mat.make()
		materials[mat.title] = mat

	# Cleanup
	del matini

#-------------------------------------------------------------------------------
# Create a new ini file only with the modified materials
#-------------------------------------------------------------------------------
def makeIni(modified=True):
	global materials, groups

	matini = configparser.RawConfigParser()

	# Create groups
	matini.add_section("groups")
	for i in range(len(groups)):
		matini.set("groups", "group.%d"%(i+1), groups[i])

	# Add materials
	for mat in list(materials.values()):
		if modified and mat.modified == "": continue

		section = mat.title
		matini.add_section(section)
#		matini.set(section, "name",     mat.name.encode("utf-8"))
		matini.set(section, "name",     mat.name)
		matini.set(section, "groups",   " ".join(map(str,mat.group)))
		matini.set(section, "composition", _COMPOSITION[mat.composition])
		matini.set(section, "elements", len(mat.elements))
		matini.set(section, "state",    _STATE[mat.state])
		if mat.modified != "":
			matini.set(section, "modified", mat.modified)

		matini.set(section, "density",  str(mat.density))
		matini.set(section, "pressure", str(mat.pressure))
		matini.set(section, "rhor",     str(mat.rhor))
		matini.set(section, "ionization", str(mat.ionization))
		matini.set(section, "dpa",      str(mat.dpa_threshold))
		matini.set(section, "xrho_dEdx",str(mat.xrho_dEdx))
		matini.set(section, "xrho",     str(mat.xrho))

		matini.set(section, "s_cbar",   str(mat.s_cbar))
		matini.set(section, "s_x0",     str(mat.s_x0))
		matini.set(section, "s_x1",     str(mat.s_x1))
		matini.set(section, "s_a",      str(mat.s_a))
		matini.set(section, "s_m",      str(mat.s_m))
		matini.set(section, "s_d0",     str(mat.s_d0))
		i = 1
		for Z,A,f in mat.elements:
			matini.set(section, "z.%d"%(i), Z)
			matini.set(section, "a.%d"%(i), A)
			matini.set(section, "f.%d"%(i), f)
			i += 1
		matini.set(section, "notes", \
			"$"+mat.notes.replace("\n","\n$"))
#		matini.set(section, "notes", \
#			"$"+string.replace(mat.notes.encode("utf-8"),"\n","\n$"))
	return matini

#===============================================================================
if __name__ == "__main__":
	import os

	rootdir = os.path.dirname(sys.argv[0])
	Input.init()
	PeriodicTable.init( os.path.join(rootdir, "db/isotopes.ini"))
	init(               os.path.join(rootdir, "db/material.ini"))

	root = Tk()
	root.deiconify()
	tkFlair.loadIcons()

	mat = MaterialFrame(root)
	mat.showFrame()
	mat.refresh()
	mat.pack(expand=YES, fill=BOTH)

	root.mainloop()
	tkFlair.delIcons()

	matini = makeIni()
	fini = open("mat.ini","w")
	matini.write(fini)
	fini.close()
