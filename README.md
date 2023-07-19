# Java Implementation of Recommendation ITU-R P.452

This code repository is a development branch of a Java software implementation of  [Recommendation ITU-R P.452-18](https://www.itu.int/rec/R-REC-P.452/en)  with a prediction procedure for the evaluation of interference between stations on the surface of the Earth at frequencies above about 0.1 GHz. 

[//]: < This code repository contains a Java software implementation of  [Recommendation ITU-R P.452-18](https://www.itu.int/rec/R-REC-P.452/en)  with a prediction procedure for the evaluation of interference between stations on the surface of the Earth at frequencies above about 0.1 GHz. >

[//]: < This version of the code is functionally identical to the reference version approved by ITU-R Working Party 3M and published by Study Group 3 on [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx). This version of the code is also implemented in [SEAMCAT](https://seamcat.org).>


The following table describes the structure of the folder `./src/` containing the Java implementation of Recommendation ITU-R P.452.

| File/Folder               | Description                                                         |
|----------------------------|---------------------------------------------------------------------|
|`main/P452.java`                | Java class implementing Recommendation ITU-R P.452-18         |
|`test/P452Test.java`          | Java class implementing validation tests against the reference MATLAB/Octave implementation of this Recommendation for a range of input variables.          |



## Function Call

~~~ 
Lb = tl_p452(maps, f, p, d, h, g, zone, htg, hrg, phit_e, phit_n, phir_e, phir_n, Gt, Gr, pol, dct, dcr, press, temp)
~~~

## Required input arguments of function `tl_p452`

| Variable          | Type   | Units | Limits       | Description  |
|-------------------|--------|-------|--------------|--------------|
| `maps`           | class `P452DigitalMaps` | |  | Object containing all the digital maps (DN50, N050) necessary for computation |
| `f`               | scalar double | GHz   | ~0.1 ≤ `f` ≤ ~50 | Frequency   |
| `p         `      | scalar double | %     | 0.001 ≤ `p` ≤ 50 | Time percentage for which the calculated basic transmission loss is not exceeded |
| `d`               | array double | km    |  0 < `max(d)` ≤ ~10000 | Terrain profile distances (in the ascending order from the transmitter)|
| `h`          | array double | m (asl)   |   | Terrain profile heights |
| `g`          | array double | m (asl)   |  | Clutter + Terrain profile heights   |
| `zone`           | array int    |       | 1 - Coastal Land, 2 - Inland, 3 - Sea             |  Radio-climatic zone types |
| `htg`           | scalar double    | m      |           |  Tx antenna height above ground level |
| `hrg`           | scalar double    | m      |          |  Rx antenna height above ground level |
| `phit_e`           | scalar double    | deg      |     0 ≤ `phit_e`  ≤ 360  or -180 ≤ `phit_e`  ≤ 180  |  Tx longitude |
| `phit_n`           | scalar double    | deg      |     -90 ≤ `phit_n`  ≤ 90          |  Tx latitude |
| `phir_e`           | scalar double    | deg      |     0 ≤ `phir_e`  ≤ 360   or  -180 ≤ `phir_e`  ≤ 180       |  Rx longitude |
| `phir_n`           | scalar double    | deg      |     -90 ≤ `phir_n`  ≤ 90          |  Rx latitude |
| `Gt`,  `Gr`           | scalar double  |   dBi    |           |  Tx/Rx antenna gain in the direction of the horizon towards along the great-circle interference path. |
| `pol`           | scalar int    |       |   `pol`  = 1, 2          |  Polarization of the signal: 1 - horizontal, 2 - vertical |
| `dct`           | scalar double    | km      |   `dct` ≥ 0          |  Distance over land from the Tx antenna to the coast along the great-circle interference path. To be set to zero for a terminal on a ship or sea platform.|
| `dcr`           | scalar double    | km      |   `dcr` ≥ 0          |  Distance over land from the Rx antenna to the coast along the great-circle interference path. To be set to zero for a terminal on a ship or sea platform.|
| `press`           | scalar double    | hPa      |             | Dry air pressure.|
| `temp`           | scalar double    | deg C      |             | Air temperature.|


 
## Outputs ##

| Variable   | Type   | Units | Description |
|------------|--------|-------|-------------|
| `Lb`    | double | dB    | Basic transmission loss |




## References

* [Recommendation ITU-R P.452](https://www.itu.int/rec/R-REC-P.452/en)

* [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx)

* [MATLAB/Octave Implementation of Recommendation ITU-R P.452](https://github/eeveetza/p452)

* [SEAMCAT - Spectrum Engineering Advanced Monte Carlo Analysis Tool](https://seamcat.org)