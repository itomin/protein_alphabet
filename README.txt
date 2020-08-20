ssh root@widdas01 -i M:\.ssh\id_rsa -L 8888:localhost:8888
ssh ubuntu@3.121.229.250 -p 8080  -i biot.pem -v -L 8888:localhost:8888

sudo docker run --rm -it -p 8888:8888 -v /mnt:/home/nb ocramz/jupyter-docker-pymol
docker run --rm -it -p 8888:8888 -v dev:/home/nb pymol:pandas



pymol_render("1xxm", 1013, 1027, "1013-1027/")
pymol_render("3gmv", 13, 27, "2-157/")

fetch 1xxm
color palegreen
fetch 3gmv
color lightblue
select frag1_sele, 1xxm and resi 1013-1027
select frag2_sele, 3gmv and resi 13-27

color green, frag1_sele
color blue, frag2_sele
super frag1_sele, frag2_sele
copy_to frag1_sele, frag1
copy_to frag2_sele, frag2
hide everything
show cartoon
