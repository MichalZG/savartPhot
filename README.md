# savartPhot

!!! ONLY FOR 64-bit !!!
1. conda create --name savart --file spec-file.txt
2. Create link to program:
  * sudo ln -s path_to_run_script /usr/bin/name
  * sudo chmod +x /usr/bin/name
3. Start program with typing name 


Examples:

* run list 
python run_list_objects.py /home/pi/Programs/python-programs/savartPhot/var/pipeline_config_files/config_list_objects.cfg /home/pi/savarts/dane_suhora/2017-11-16/all/
* run stack
python run_stack.py /home/pi/Programs/python-programs/savartPhot/var/pipeline_config_files/config_stack.cfg /home/pi/savarts/dane_suhora/2017-11-16/hd204827/ 8 /home/pi/savarts/dane_suhora/2017-11-16/hd204827/output

* photometry need coordinate file of stars, where y1 > y2 !

P1 x1 y1 x2 y2

P3 x1 y1 x2 y2\\

example:

P1 544 233 450 245

P3 566 344 477 344

* run photometry
python run_savart_phot.py /home/pi/Programs/python-programs/savartPhot/var/pipeline_config_files/config_savartphot.cfg /home/pi/savarts/dane_suhora/2017-11-16/hd204827/output/ /home/pi/savarts/dane_suhora/2017-11-16/hd204827/coo.dat /home/pi/savarts/dane_suhora/2017-11-16/hd204827/phot_output/


