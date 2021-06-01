videos: hpprp
	make -B density/density.avi >/dev/null
	make -B rest/rest.avi >/dev/null
	make -B move/move.avi >/dev/null

hpprp: ucodes.cpp collide.inc
	luna --verbose --no-cleanup test.fa 2>err


density/density.avi:
	system/image_gen.sh density density/density.avi

rest/rest.avi:
	system/image_gen.sh rest rest/rest.avi

move/move.avi:
	system/image_gen.sh move move/move.avi

