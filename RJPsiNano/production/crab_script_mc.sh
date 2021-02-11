echo "================= CMSRUN starting ===================="
cmsRun -j FrameworkJobReport.xml -p PSet.py 
#cmsRun -j FrameworkJobReport.xml -p ../test/run_nano_jpsi_cfg.py
echo "================= CMSRUN finished ===================="
#right name for the input file
echo "==================puReweight starting ================"
python  puReweight.py

# output file slightly different name
#rm input file
echo "================= puReweight finished ================"


