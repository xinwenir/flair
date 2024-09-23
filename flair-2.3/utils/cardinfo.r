f = open("/soft/fluka/fluka2006.manual","r")
notInCards = 1
cards. = ""
do forever
	line = read(f)
	if eof(f) then leave
	if left(line,1)=="*" then iterate
	if notInCards then do
		if line="Introduction to the FLUKA input options" then do
			do until card=="ASSIGNMAt"
				line = read(f)
				parse var line card text
			end
			parse var line card text
			fullname = strip(translate(card))
			card = strip(left(card,8))
			text = strip(text)
			text = translate(left(text,1))||substr(text,2)
			cards.card = strip(text)
			cards.card.@FULL = fullname
			notInCards = 0
		end
	end; else do
		if line = "" then leave
		if left(line,11)="" then
			cards.card = cards.card strip(line)
		else do
			parse var line card text
			fullname = strip(translate(card))
			card = strip(left(card,8))
			text = strip(text)
			cards.card = translate(left(text,1))||substr(text,2)
			cards.card.@FULL = fullname
		end
	end
end
call close f

f     = open("cardorder.csv","r")
fout  = open("__cardinfo","w")
order = 0
do forever
	line = read(f)
	if eof(f) then leave
	if left(line,1)=="#" | line="" then iterate
	parse var line tag categories desc
	parse var tag '"' tag '"'
	categories = '("' || changestr(', ',categories,'", "') || '")'
	say tag
	if cards.tag.@FULL="" then cards.tag.@FULL = tag
	call lineout fout,'_addInfo("'cards.tag.@FULL'",'
	call lineout fout,'	'order","
	call lineout fout,'	'categories','
	if desc=="" then
		call lineout fout,'	"'changestr('"',cards.tag,"'")'",'
	else
		call lineout fout,'	"'changestr('"',desc,"'")'",'
	call lineout fout,'	_'changestr("-",tag,"_")'_layout)'
	order = order + 1
end
call close f
call close fout

say "__cardinfo.py file is generated"
