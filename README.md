# savartPhot

!!! ONLY FOR 64-bit !!!
1. conda create --name savart --file spec-file.txt
2. Create link to program:
  * sudo ln -s path_to_run_script /usr/bin/name
  * sudo chmod +x /usr/bin/name
3. Start program with typing name 


Examples:

python run_list_objects.py /home/pi/Programs/python-programs/savartPhot/var/pipeline_config_files/config_list_objects.cfg /home/pi/savarts/dane_suhora/2017-11-16/all/

python run_stack.py /home/pi/Programs/python-programs/savartPhot/var/pipeline_config_files/config_stack.cfg /home/pi/savarts/dane_suhora/2017-11-16/hd204827/ 8 /home/pi/savarts/dane_suhora/2017-11-16/hd204827/output

python run_savart_phot.py /home/pi/Programs/python-programs/savartPhot/var/pipeline_config_files/config_savartphot.cfg /home/pi/savarts/dane_suhora/2017-11-16/hd204827/output/ /home/pi/savarts/dane_suhora/2017-11-16/hd204827/coo.dat /home/pi/savarts/dane_suhora/2017-11-16/hd204827/phot_output/
