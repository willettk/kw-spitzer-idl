pro kwphot, fcb_0, fcb_1, pos1, hdr_1

; NOTE: mysmart must already be invoked from the command line for this program to run

;filedir = '~/Astronomy/Research/Spitzer/ohm/data/mega0??_peakup/ch0/bcd/*bcd*fits'
filedir = '~/Astronomy/Research/Spitzer/ohm/data/mega002_peakup/ch0/bcd/*bcdb*fits'
files = file_search(filedir)

fcb_0 = readfits(files(0),hdr_0)
fcb_1 = readfits(files(1),hdr_1)
fcb_2 = readfits(files(2),hdr_2)
fcb_3 = readfits(files(3),hdr_3)
fcb_4 = readfits(files(4),hdr_4)
fcb_5 = readfits(files(5),hdr_5)

pos1 = fcb_0 < fcb_1
pos2 = fcb_2 < fcb_3
pos3 = fcb_4 < fcb_5

;stop
end
