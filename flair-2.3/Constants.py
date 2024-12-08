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
# Date:	27-Jan-2012

# Bit flags for object selection
# XXX WARNING should be the same as in the geoviewer XXX

BIT_SELECT	= 1
BIT_LOCK	= 1<<1
BIT_FREEZE	= 1<<2
BIT_VISIBLE	= 1<<3
BIT_WIREFRAME   = 1<<4
BIT_MOVE        = 1<<5
BIT_BBOX        = 1<<6

BIT_SELFREEZE   = BIT_SELECT | BIT_FREEZE

BIT_UNSELECT    = ~BIT_SELECT
BIT_UNLOCK      = ~BIT_LOCK
BIT_UNFREEZE    = ~BIT_FREEZE
BIT_UNVISIBLE   = ~BIT_VISIBLE
BIT_UNWIREFRAME = ~BIT_WIREFRAME
BIT_UNBBOX      = ~BIT_BBOX

# XXX WARNING should be the same as in the geoviewer XXX
DRAW_CLEAR     = 1
DRAW_SEGMENTS  = 1<<1
DRAW_REGIONS   = 1<<2
DRAW_GRID      = 1<<3
DRAW_AXES      = 1<<4
DRAW_TITLE     = 1<<5
DRAW_STATUS    = 1<<6
DRAW_LABELS    = 1<<7
DRAW_IMAGE     = 1<<8
DRAW_USRBIN    = 1<<9
DRAW_3D        = 1<<10
DRAW_EDGE      = 1<<11
DRAW_VERTICES  = 1<<12
DRAW_MGDRAW    = 1<<13
DRAW_COLORBOX  = 1<<14
DRAW_LATTICES  = 1<<15
DRAW_VOXEL     = 1<<16
DRAW_WIREFRAME = 1<<17
DRAW_BBOX      = 1<<18
DRAW_MESSAGE   = 1<<19
DRAW_TEST      = 1<<30

DRAW_FAST      = DRAW_3D|DRAW_WIREFRAME|DRAW_BBOX|DRAW_CLEAR

ALPHA    = "alpha"
ID       = "@id"
SELECT   = "@select"
COLOR    = "@color"
MATERIAL = "@material"

VIEWNAME   = ["Red",   "Green",  "Blue",   "Magenta"]
VIEWCOLOR  = [0xFF7070, 0x70FF70, 0x7070FF, 0xFF70FF]

MAXERRORS  = 50
MAXHISTORY = 20			# Save last views

# Actions
# Negative marks the sub-task, that do not hide the task dialog

# Selections
ACTION_SELECT        =   0	# Default action (select)
ACTION_AREA_SELECT   =   1	# Select bodies or regions using rectangle
ACTION_SEL_OR_MOVE   =   2	# Still fuzzy if we want to select or to move
ACTION_SEL_OR_AREA   =   3	# Still fuzzy if we want to select on point or an area
ACTION_SELECT_REGION =  -4	# Select a region
ACTION_SELECT_ZONE   =  -5	# Select an existing zone from a region
ACTION_SELECT_BODIES =  -6	# Select bodies

# Display
ACTION_PAN           =  10	# Pan display
ACTION_PAN_SLOWDOWN  =  11	# Pan slowdown
ACTION_ORBIT         =  12	# Orbit (trackball) movement
ACTION_ZOOM_IN       =  13
ACTION_ZOOM_OUT      =  14
ACTION_ZOOM_ON       =  15

# Viewport lines
ACTION_VIEW_CENTER   =  20	# Move viewport center
ACTION_VIEW_MOVE     =  21	# Move viewport line leaving center unchanged in other axis
ACTION_VIEW_ROTATE   =  22	# Rotate viewport line
ACTION_VIEW_CENTER   =  23       # Center viewports to location

# Addition
ACTION_ADD           =  30	# Add an object
ACTION_ADD_NEXT      = -31	# Add next point of object
ACTION_ADD_ZONE      = -32
ACTION_EDIT          =  33	# Edit/Correct overlaps
# Objects manipulation
ACTION_MOVE          =  34	# Move selected objects
ACTION_ROTATE        =  35	# Rotate selected objects

# Other editing
ACTION_ZONE          =  40	# Select zone master action
ACTION_ZONEPAINT     =  42	# Select multiple zones
ACTION_PEN           =  43	# Draw using pen
ACTION_PAINT         =  44	# Paint region properties
ACTION_INFO          =  45	# Display info at cursor
