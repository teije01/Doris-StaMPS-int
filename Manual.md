## Manual for Doris StaMPS interface
Short description of how to use these scripts. 
# Installation
The following programs need to be installed:
* [Dorisv5](http://doris.tudelft.nl) 
* [Triangle](http://www.cs.cmu.edu/˜quake/triangle.html)
* [Snaphu](http://www.stanford.edu/group/radar/softwareandlinks/sw/snaphu/)
* [StaMPS](https://homepages.see.leeds.ac.uk/~earahoo/stamps/) v3.3b1

StaMPS is installed with the following commands:
```shell
tar -xvf StaMPS_v3.3b1.tar
cd StaMPS_v3.3b1/src
make
make install
```

In the StaMPS_v3.3b1/ folder, extract the Sentinel add-on and make it executable by entering:
```shell
cd StaMPS_v3.3b1/Sentinel/
chmod 750 *
```

Now edit the StaMPS_CONFIG.bash and include at least the following:
(replace user by actual user)
```bash
#NECESSARY
export SOFT="/home/user/STAMPS"
export STAMPS="$SOFT/StaMPS_v3.3b1"
export TRIANGLE_BIN="$SOFT/triangle"
export SNAPHU_BIN="/home/user/bin/snaphu"
export DORIS_BIN="/home/user/bin/doris/doris_current"
export S1ST_DIR="$STAMPS/Sentinel"

#####################################
# shouldn't need to change below here
#####################################
export MATLABPATH=$STAMPS/matlab:`echo $MATLABPATH`
export DORIS_SCR="$STAMPS/DORIS_SCR"

#The command below is on one line and needs to contain all folders that need to be added to the search path (adduming doris is already added to the search path)
export PATH=${PATH}:$STAMPS/bin:$DORIS_SCR:$TRIANGLE_BIN:$SNAPHU_BIN:$S1ST_DIR
```

To activate its contents, execute the following line:
```shell
source StaMPS_CONFIG.bash
```
and consider adding it to your .bashrc so it is activated upon login
```shell
source /home/path_to_StaMPS_v3.3b1/StaMPS_CONFIG.bash
```

# Script overview (short)
After creating a single master-stack with Doris, the interface can be used to get StaMPS up and running using the following scripts:
1. `link_ifgs_sentinel_master`
   calls to:
   * `get_slc_lp_S1`.
2. `link_res_sentinel_slaves`
3. `crop_master_slc_dem`
   calls to:
   * `CropData.py` 
   * `edit_res_crop_master`
4. crop_slave_slc
   calls to:
   * `CropData.py` 
   * `edit_res_crop_slave`
5. `mt_prep_sentinel`
  calls to:
   * `mt_extract_info_sentinel`
   * `mt_bperp_angle_sentinel` 
      calls to:
      * `LookAngle.py`
      * `BperpDate.py`
   * `mt_extract_cands`

# How to use the scripts:

Work in progress
