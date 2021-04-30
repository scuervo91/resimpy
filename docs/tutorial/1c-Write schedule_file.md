# Write Schedule Data File

The `schedule_writer` is a simple function to write into a string the keywords and data that depends on the time. It includes 'DATES', 'WELSPECS', 'COMPDAT', 'GRUPTREE' or any keyword that needs to be placed in the right time.


```python
import pandas as pd
from resimpy import schedule_writer
```

For this example, we can import the keywords from a file to a pandas.DataFrame. It assumes the dataframe contains a column named 'date'. 

The proccess consist on iterate over all dates contained in all available keywords and write the keywords according the date. 


```python
all_keywords_toexport = {}
some_keywords = ['DATES', 'WELSPECS', 'COMPDAT', 'GRUPTREE']
for keyword in some_keywords:
    all_keywords_toexport[keyword] = pd.read_excel('prod_input_test.xlsx',sheet_name=keyword)

    print(keyword)
    print(all_keywords_toexport[keyword].head(3))
    print('-----------------------------')
```

    DATES
       DAY MONTH  YEAR       date
    0   31   OCT  2022 2022-10-31
    1    8   NOV  2022 2022-11-08
    2   25   NOV  2022 2022-11-25
    -----------------------------
    WELSPECS
       WELL     GROUP   I   J BHPREF TYPE DRADIUS INFLOW  AUTO XFLOW PVTNUM  \
    0  D-1H   MANI-D1  22  22     1*  OIL      1*     1*  STOP    1*     1*   
    1  B-2H  B1-DUMMY  15  31     1*  OIL      1*     1*  STOP    1*     1*   
    2  D-2H  D2-DUMMY  14  28     1*  OIL      1*     1*  STOP    1*     1*   
    
      DENOPT FIPNUM       date  
    0     1*     1* 2022-10-31  
    1     1*     1* 2022-10-31  
    2     1*     1* 2022-12-11  
    -----------------------------
    COMPDAT
       WELL   I   J  K1  K2 Status Sat.     CF   DIAM       KH SKIN  ND DIR  \
    0  D-1H  22  22   5   5   OPEN   1*  5.505  0.216  510.312   1*  1*   Z   
    1  D-1H  23  22   6   6   OPEN   1*  0.101  0.216    9.456   1*  1*   Z   
    2  D-1H  23  22   7   7   OPEN   1*  4.938  0.216  452.905   1*  1*   Z   
    
           Ro       date  
    0  15.511 2022-10-31  
    1  16.532 2022-10-31  
    2  14.704 2022-10-31  
    -----------------------------
    GRUPTREE
         LOWER HIGHER       date
    0     INJE  FIELD 2022-10-31
    1     PROD  FIELD 2022-10-31
    2  MANI-B2   PROD 2022-10-31
    -----------------------------



```python
sched_string = schedule_writer(all_keywords_toexport)

with open('schedule_test.inc','w') as file:
    file.write(sched_string)
```

Print first part of the file 

```
RPTRST
 'BASIC=1' /
RPTSCHED
 'FIP=1' 'WELLS=1' 'WELLSPECS' /
DATES
31 'OCT' 2022 /
/
WELSPECS
D-1H MANI-D1 22 22 1* OIL 1* 1* STOP 1* 1* 1* 1*/
B-2H B1-DUMMY 15 31 1* OIL 1* 1* STOP 1* 1* 1* 1*/
C-4H MANI-C 11 35 1* OIL 1* 1* 1* 1* 1* 1* 1*/
/
COMPDAT
D-1H 22 22 5 5 OPEN 1* 5.505 0.216 510.312 1* 1* Z 15.511/
D-1H 23 22 6 6 OPEN 1* 0.101 0.216 9.456 1* 1* Z 16.532/
D-1H 23 22 7 7 OPEN 1* 4.938 0.216 452.905 1* 1* Z 14.704/
D-1H 23 22 9 9 OPEN 1* 19.086 0.216 1745.284 1* 1* Z 14.493/
D-1H 23 22 10 10 OPEN 1* 50.101 0.216 4655.453 1* 1* Z 15.689/
D-1H 23 22 11 11 OPEN 1* 8.974 0.216 823.585 1* 1* Z 14.751/
D-1H 23 22 12 12 OPEN 1* 0.479 0.216 43.304 1* 1* Z 13.707/
D-1H 23 22 13 13 OPEN 1* 12.603 0.216 1152.42 1* 1* Z 14.489/
B-2H 17 31 9 9 OPEN 1* 17.246 0.216 1285.863 1* 1* X 5.865/
B-2H 19 31 9 9 OPEN 1* 13.2 0.216 991.575 1* 1* X 6.044/
B-2H 20 31 10 10 OPEN 1* 36.54 0.216 2804.161 1* 1* X 6.593/
B-2H 21 31 10 10 OPEN 1* 12.052 0.216 921.178 1* 1* X 6.486/
B-2H 22 32 10 10 OPEN 1* 67.732 0.216 5174.542 1* 1* X 6.472/
B-2H 24 32 10 10 OPEN 1* 42.421 0.216 3232.419 1* 1* X 6.404/
B-2H 25 32 10 10 OPEN 1* 29.697 0.216 2261.93 1* 1* X 6.393/
B-2H 29 33 10 10 OPEN 1* 10.49 0.216 807.533 1* 1* X 6.677/
C-4H 11 35 1 1 OPEN 1* 45.314 0.216 4253.571 1* 1* Z 16.503/
C-4H 11 35 2 2 OPEN 1* 43.674 0.216 4103.809 1* 1* Z 16.588/
/
GRUPTREE
INJE FIELD/
PROD FIELD/
MANI-B2 PROD/
MANI-B1 PROD/
MANI-D1 PROD/
MANI-D2 PROD/
MANI-E1 PROD/
MANI-E2 PROD/
MANI-K1 MANI-B1/
MANI-K2 MANI-D2/
MANI-C INJE/
MANI-F INJE/
WI-GSEG INJE/
B1-DUMMY MANI-B1/
D2-DUMMY MANI-D2/
/
DATES
08 'NOV' 2022 /
25 'NOV' 2022 /
11 'DEC' 2022 /
/
WELSPECS
D-2H D2-DUMMY 14 28 1* OIL 1* 1* STOP 1* 1* 1* 1*/
/
COMPDAT
D-2H 14 26 9 9 OPEN 1* 21.45 0.216 1590.754 1* 1* Y 5.741/
D-2H 14 25 9 9 OPEN 1* 39.557 0.216 2921.561 1* 1* Y 5.648/
D-2H 14 23 9 9 OPEN 1* 10.183 0.216 748.871 1* 1* Y 5.554/
D-2H 14 22 9 9 OPEN 1* 121.842 0.216 8821.805 1* 1* Y 5.225/
D-2H 14 21 9 9 OPEN 1* 140.551 0.216 10196.747 1* 1* Y 5.266/
D-2H 14 20 9 9 OPEN 1* 24.486 0.216 1793.318 1* 1* Y 5.465/
D-2H 14 15 9 9 OPEN 1* 29.883 0.216 2344.667 1* 1* Y 7.229/
D-2H 14 14 9 9 OPEN 1* 82.852 0.216 6372.295 1* 1* Y 6.653/
D-2H 14 13 9 9 OPEN 1* 24.664 0.216 1809.697 1* 1* Y 5.504/
/
DATES
26 'DEC' 2022 /
26 'JAN' 2023 /
23 'FEB' 2023 /
23 'MAR' 2023 /
/
WELSPECS
B-4H B1-DUMMY 10 32 1* OIL 1* 1* STOP 1* 1* 1* 1*/
/
```
