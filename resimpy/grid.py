#########################################################################
# Most of the Code has been taken from the next  Github Repository:      #
#  https://github.com/BinWang0213/PyGRDECL                              #
#  Code is used to load and manipulate Eclipse Data Grid
#########################################################################

import numpy as np 
import pyvista as pv 
import vtk
from shapely.geometry import Point
import math
import os
import pandas as pd 
from pydantic import BaseModel, Field, Extra, validator
from enum import Enum
from typing import List, Optional, Union, Dict

petrophysical_properties = ['PORO','PERMX','PERMY','PERMZ','SW','RT']

SupportKeyWords=[
    'SPECGRID', #Dimenion of the corner point grid
    'DIMENS',   #Define the dimension of the cartesian grid
    'TOPS','DX','DY','DZ',
    'COORD','ZCORN',
    'PORO',
    'PERMX' , 'PERMXY', 'PERMXZ', 
    'PERMYX', 'PERMY' , 'PERMYZ', 
    'PERMZX', 'PERMZY', 'PERMZ',
    'ACTNUM',
    'SATNUM', 'NTG',
    'INCLUDE',
    
]

KeyWordsDatatypes=[#Corrsponding data types
    int,
    int,
    int,int,int,int,
    float,float,
    float,
    float,float,float,
    float,float,float,
    float,float,float,
    int,
    int,float
]

def parseDataArray(DataArray):
        """Parse special dataArray format in GRDECL 
        example:
            5*3.0=[3.0 3.0 3.0 3.0 3.0]
            1.0 2*3.0 5.0=[1.0 3.0 3.0 5.0]
        
        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """

        data=[]
        error_count=0
        for value in DataArray:
            if(is_number(value)==2):
                num,val=value.split('*')
                for i in range(int(num)): data.append(val)
            elif(is_number(value)==1):
                data.append(value)
            else:
                error_count+=1
        
        if(error_count>0):
            print(DataArray)
        
        assert error_count==0, '[Error] Can not find any numeric value!'
        
        return data
    
def is_number(s):
    #Determine a string is a number or not
    #Used in [read_GRDECL] [getBlkdata]
    try:
        float(s)
        return True
    except ValueError:
        pass
 
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    
    try: #Special format N*val= [val val val ....]
        num, val = s.split('*')
        return 2
    except ValueError:
        pass
 
    return False

## Auxilary functions
def cell_id(i,j,k,nx,ny):
    """
    Get the cell Id given i,j,k indexes. 
        * ---  *  ---  *  --- *
        | 0,2  |  1,2  |  2,2 |  <- Cell id 6,7,8
        * ---  *  ---  *  --- *
        | 0,1  |  1,1  |  2,1 |  <- Cell id 3,4,5
        * ---  *  ---  *  --- *
        | 0,0  |  1,0  |  2,0 |  <- Cell id 0,1,2
        * ---  *  ---  * ---  *
    """
    cell = (nx*j+i)+k*nx*ny

    return cell

def cell_ijk(cell_id,nx,ny):
    """
    Get the cell indexes i,j,k given the cell id
        * ---  *  ---  *  --- *
        | 0,2  |  1,2  |  2,2 |  <- Cell id 6,7,8
        * ---  *  ---  *  --- *
        | 0,1  |  1,1  |  2,1 |  <- Cell id 3,4,5
        * ---  *  ---  *  --- *
        | 0,0  |  1,0  |  2,0 |  <- Cell id 0,1,2
        * ---  *  ---  * ---  *
    """
    if cell_id==0:
        return 0,0,0
    k=math.ceil(cell_id/(nx*ny))-1
    j=math.ceil((cell_id-(nx*ny)*k)/nx)-1
    i=math.ceil(cell_id-(nx*ny*k)-nx*j)
    return i,j,k

#Interpolate z on pillars
def interpolate_z_pillar(z,p):
    """
    Obtain the eight coords for a cell
        X,Y coords has to be interpolated from Z
    xy1=xy0+k*z
    Pillar=np.array([[x0 y0 z0],[x1 y1 z1]])
    """
    x = ((p[0,0]-p[1,0])/(p[0,2]-p[1,2]))*(z-p[0,2])+p[0,0]
    y = ((p[0,1]-p[1,1])/(p[0,2]-p[1,2]))*(z-p[0,2])+p[0,1]

    xyz = np.array([x,y,z])
    return xyz

#3D Rotation funtions
def rotation(points,azimuth,dip,plunge):
    assert points.ndim == 2

    azi_rad = np.radians(azimuth)
    dip_rad = np.radians(dip)
    plg_rad = np.radians(plunge)

    ry = np.array([
        [np.cos(plg_rad),0,-np.sin(plg_rad)],
        [0,1,0],
        [np.sin(plg_rad),0,np.cos(plg_rad)],
    ])
    rx = np.array([
        [1,0,0],
        [0,np.cos(dip_rad),np.sin(dip_rad)],
        [0,-np.sin(dip_rad),np.cos(dip_rad)]
    ])
    rz = np.array([
        [np.cos(azi_rad),-np.sin(azi_rad),0],
        [np.sin(azi_rad),np.cos(azi_rad),0],
        [0,0,1]
    ])

    rot = np.matmul(np.matmul(ry,rx),rz)

    rot_points = np.matmul(points,rot)

    return rot_points

def RemoveCommentLines(data,commenter='--'):
    #Remove comment and empty lines
    data_lines=data.strip().split('\n')
    newdata=[]
    for line in data_lines:
        if line.startswith(commenter) or not line.strip():
            # skip comments and blank lines
            continue   
        newdata.append(line)
    return '\n'.join(newdata)


def scanKeyword(data):
    #scan and find the keyword
    #e.g. ['INIT','DX','2500*30.0'] -> ['DX','2500*30.0']
    for key in SupportKeyWords:
        if (key in data) and (data.find(key)!=0):
            return data[data.find(key):-1]
    return data
## Grid Class
class GridTypeEnum(str, Enum): 
    cartesian = 'cartesian'
    corner_point = 'corner_point'
    
class Grid(BaseModel):
    """
    Class for Reservoir Simulation Grid 
        * Cartesian
        * Corner-Point
    """
    grid_type: GridTypeEnum = Field(None)
    nx: int = Field(None, gt=0)
    ny: int = Field(None, gt=0)
    nz: int = Field(None, gt=0)
    tops: Optional[Union[List[float],float]] = Field(None)
    dx: Optional[Union[List[float],float]] = Field(None)
    dy: Optional[Union[List[float],float]] = Field(None)
    dz: Optional[Union[List[float],float]] = Field(None)
    origin: Point = Field(Point(0,0))
    azimuth: float = Field(0, ge=0, le=360)
    dip: float = Field(0, ge=0, le=90)
    plunge: float = Field(0, ge=0, le=90)
    coord: Optional[List[float]] = Field(None)
    zcorn: Optional[List[float]] = Field(None)
    spatial_data: Optional[Dict[str,Union[List[float],float]]] = Field(None)
    skiped_keywords: int = Field(0)
    class Config:
        arbitrary_types_allowed = True
        validate_assignment = True
        extra = Extra.forbid
        
    @validator('tops')
    def match_tops(cls,v,values):
        if values['grid_type'] != GridTypeEnum.cartesian:
            raise ValueError('Deltas must be set on cartesian grid')
        n = values['nx'] * values['ny'] * values['nz']
        ntop = values['nx'] * values['ny']
        if isinstance(v,list):
            assert any([len(v) == n,  len(v) == ntop])
            return v 
        else:     
            return np.full(ntop,v).tolist()
        
    @validator('dx','dy','dz')
    def match_deltas(cls,v,values):
        if values['grid_type'] != GridTypeEnum.cartesian:
            raise ValueError('Deltas must be set on cartesian grid')
        n = values['nx'] * values['ny'] * values['nz']
        if isinstance(v,list):
            assert len(v) == n 
            return v 
        else:     
            return np.full(n,v).tolist()
        
   
    @validator('origin')
    def check_origin_coord(cls,v,values):
        if values['grid_type'] != GridTypeEnum.cartesian:
            raise ValueError('Origin must be set only on castesian grid')
        assert v.has_z
        return v
    
    
    @validator('coord')
    def check_coord(cls,v,values):
        if values['grid_type'] == GridTypeEnum.cartesian:
            raise ValueError('Coord must be set on corner point grid')
        length = 6*(values['nx']+1)*(values['ny']+1)
        assert len(v) == length, f'list must be of length{length}'
        return v
    
    @validator('zcorn')
    def check_zcorn(cls,v,values):
        if values['grid_type'] == GridTypeEnum.cartesian:
            raise ValueError('Coord must be set on corner point grid')
        length = 8 * values['nx'] * values['ny'] * values['nz']
        assert len(v) == length, f'list must be of length{length}'
        return v
    
    @validator('spatial_data')
    def check_spatial_data(cls,v,values):
        n = values['nx'] * values['ny'] * values['nz']
        for i in v:
            if isinstance(v[i],list):
                assert len(v[i]) == n
                return v 
            else:     
                return np.full(n,v[i]).tolist()
    
    @property
    def n(self):
        return self.nx * self.ny * self.nz

    def cartesian_vertices_coord(self):
        #Vertices coordinates starting at 0,0,0
        x_vert_cord = np.concatenate((np.zeros(1),np.array(self.dx).reshape((self.nx,self.ny,self.nz),order='f')[:,0,0]),axis=0).cumsum()
        y_vert_cord = np.concatenate((np.zeros(1),np.array(self.dy).reshape((self.nx,self.ny,self.nz),order='f')[0,:,0]),axis=0).cumsum()
        z_vert_cord = -np.concatenate((np.zeros(1),np.array(self.dz).reshape((self.nx,self.ny,self.nz),order='f')[0,0,:]),axis=0).cumsum()

        points = np.zeros(((self.nx+1)*(self.ny+1)*(self.nz+1),3))
        for k in range(self.nz+1):
            for j in range(self.ny+1):
                for i in range(self.nx+1):
                    l = cell_id(i,j,k,self.nx+1,self.ny+1)
                    points[l,0] = x_vert_cord[i]
                    points[l,1] = y_vert_cord[j]
                    points[l,2] = z_vert_cord[k]

        #Get rotated points with respect 0,0,0
        rot_points = rotation(points,self.azimuth,self.dip,self.plunge)
        
        #Adjust the coordinates according with Origin Point
        origin = np.array([self.origin.x,self.origin.y,self.origin.z])

        _vertices_coord = rot_points + origin
        return _vertices_coord

    def cartesian_center_point_coord(self):
        #Vertices coordinates starting at 0,0,0
        x_vert_cord = np.concatenate((np.zeros(1),np.array(self.dx).reshape((self.nx,self.ny,self.nz),order='f')[:,0,0]),axis=0).cumsum()
        y_vert_cord = np.concatenate((np.zeros(1),np.array(self.dy).reshape((self.nx,self.ny,self.nz),order='f')[0,:,0]),axis=0).cumsum()
        z_vert_cord = -np.concatenate((np.zeros(1),np.array(self.dz).reshape((self.nx,self.ny,self.nz),order='f')[0,0,:]),axis=0).cumsum()

        center = np.zeros(((self.nx)*(self.ny)*(self.nz),3))
        for k in range(self.nz):
            for j in range(self.ny):
                for i in range(self.nx):
                    l = cell_id(i,j,k,self.nx,self.ny)
                    center[l,0] = np.mean((x_vert_cord[i], x_vert_cord[i+1]))
                    center[l,1] = np.mean((y_vert_cord[j], y_vert_cord[j+1]))
                    center[l,2] = np.mean((z_vert_cord[k], z_vert_cord[k+1]))

        #Get rotated points with respect 0,0,0
        rot_points = rotation(center,self.azimuth,self.dip,self.plunge)
        
        #Adjust the coordinates according with Origin Point
        origin = np.array([self.origin.x,self.origin.y,self.origin.z])

        _center_coord = rot_points + origin
        return _center_coord


#####################################################
############## Methods ###########################

    def add_spatial_data(self,key,array):
        array = np.atleast_1d(array).flatten(order='F')
        assert self.n == array.shape[0]
        spatial_dict = {key:array.tolist()}
        if self.spatial_data is None:
            self.spatial_data = spatial_dict
        else:
            self.spatial_data.update(spatial_dict)

    def read_IncludeFile(self,filename_include,NumData):
        """Read Include data file
        this data file just a series of values
        e.g. 0.2 0.3 12.23 ....
        
        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Aug. 2018
        """

        f=open(filename_include)
        contents=f.read()
        block_dataset=contents.strip().split() #Sepeart input file by slash /
        block_dataset=np.array(block_dataset,dtype=float)
        if(len(block_dataset)!=NumData):
            print('Data size %s is not equal to defined block dimension (NX*NY*NZ) %s'%(len(block_dataset),NumData))
        return block_dataset
    
    def LoadVar(self,Keyword,DataArray,DataSize):
        """Load varables into class
        example:
        
        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        if(Keyword in SupportKeyWords):#KeyWords Check
            assert len(DataArray)==DataSize,'\n     [Error-%s] Incompatible data size! %d-%d' %(Keyword,len(DataArray),DataSize)
            KeywordID=SupportKeyWords.index(Keyword)
            print('     [%s] '%(Keyword),end='')
            self.add_spatial_data(Keyword,np.array(DataArray,dtype=KeyWordsDatatypes[KeywordID]))
        else:
            print('     [Warnning] Unsupport keywords[%s]' % (Keyword))
            self.skiped_keywords+=1
    
    def read_GRDECL(self,file):
        """Read input file(GRDECL) of Reservoir Simulator- Petrel (Eclipse)  
        file format:http://petrofaq.org/wiki/Eclipse_Input_Data
        
        Arguments
        ---------
        NX, NY, NZ -- Grid dimension.
        blockData_raw -- [0] Keywords [1] values
        
        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2017
        """
        debug=0

        print('[Input] Reading ECLIPSE/PETREL file \"%s\" ....'%(file))

        #Read whole file into list
        f=open(file)
        contents=f.read()
        contents=RemoveCommentLines(contents,commenter='--')
        contents_in_block=contents.strip().split('/') #Sepeart input file by slash /
        contents_in_block = [x for x in contents_in_block if x]#Remove empty block at the end
        NumKeywords=len(contents_in_block)
        print(f'Num Keywords {NumKeywords}')
        GoodFlag=0
        for i,block in enumerate(contents_in_block):#Keyword, Block-wise
            #Clean the data where no spliter \ provided
            block=scanKeyword(block)

            blockData_raw=block.strip().split()
            Keyword=''
            DataArray=[]
            if(len(blockData_raw)>1):
                if(blockData_raw[0]=='ECHO'): #This keyword may next to real keyword
                    Keyword,DataArray=blockData_raw[1],blockData_raw[2:]                    
                else:
                    Keyword,DataArray=blockData_raw[0],blockData_raw[1:]

            #Read Grid Dimension [SPECGRID] or [DIMENS] 
            print(Keyword)
            if(Keyword=='DIMENS'):
                DataArray=np.array(DataArray[:3],dtype=int)
                self.grid_type='cartesian'
                self.nx,self.ny,self.nz=DataArray[0],DataArray[1],DataArray[2]
                print("     Grid Type=%s Grid" %(self.grid_type))
                print("     Grid Dimension(NX,NY,NZ): (%s x %s x %s)"%(self.nx,self.ny,self.nz))
                print("     NumOfGrids=%s"%(self.n))
                print('     NumOfKeywords=%s'%(NumKeywords))
                print("     Reading Keyword %d [%s] " %(i+1,Keyword),end='')
                GoodFlag=1
                continue
            elif(Keyword=='SPECGRID'):
                DataArray=np.array(DataArray[:3],dtype=int)
                self.grid_type='corner_point'
                self.nx,self.ny,self.nz=DataArray[0],DataArray[1],DataArray[2]
                print("     Grid Type=%s" %(self.grid_type))
                print("     Grid Dimension(NX,NY,NZ): (%s x %s x %s)"%(self.nx,self.ny,self.nz))
                print("     NumOfGrids=%s"%(self.n))
                print('     NumOfKeywords=%s'%(NumKeywords))
                print("     Reading Keywords [%s] " %(Keyword),end='')
                GoodFlag=1
                continue
            
            if(self.grid_type is None):#Skip unnecessary keywords
                continue

            if(Keyword in SupportKeyWords): #We need parse the special format in 
                if Keyword == 'INCLUDE':
                #if(len(DataArray)==1 and '.' in DataArray[0]):
                    folder_name=os.path.dirname(file)
                    self.read_GRDECL(os.path.join(folder_name,DataArray[0].replace("'","")))
                    continue
                    #DataArray=self.read_IncludeFile(os.path.join(folder_name,DataArray[0]),self.n)
                print(f'------{Keyword}------')

                DataArray=parseDataArray(DataArray)
            

                #Read Grid spatial information, x,y,z ordering
                #Corner point cell
                if(Keyword=='COORD'):# Pillar coords
                    assert len(DataArray)==6*(self.nx+1)*(self.ny+1),'[Error] Incompatible COORD data size!'
                    self.coord=np.array(DataArray,dtype=float).tolist()    
                elif(Keyword=='ZCORN'):# Depth coords
                    assert len(DataArray)==8*self.n, '[Error] Incompatible ZCORN data size!'
                    self.zcorn=np.array(DataArray,dtype=float)
                
                #Cartesian cell
                elif(Keyword=='DX'):# Grid size in X dir
                    assert len(DataArray)==self.n, '[Error] Incompatible DX data size!'
                    self.dx=np.array(DataArray,dtype=float).tolist() 
                elif(Keyword=='DY'):# Grid size in Y dir
                    assert len(DataArray)==self.n, '[Error] Incompatible DY data size!'
                    self.dy=np.array(DataArray,dtype=float).tolist() 
                elif(Keyword=='DZ'):# Grid size in Z dir
                    assert len(DataArray)==self.n, '[Error] Incompatible DZ data size!'
                    self.dz=np.array(DataArray,dtype=float).tolist() 
                elif(Keyword=='TOPS'):# TOP position
                    assert any([len(DataArray)==self.n,len(DataArray)==self.nx*self.ny]), '[Error] Incompatible TOPS data size!'
                    self.tops=np.array(DataArray,dtype=float).tolist() 

                #Read Grid Properties information
                else:
                    self.LoadVar(Keyword,DataArray,DataSize=self.n)

        f.close()
        #assert GoodFlag==1,'Can not find grid dimension info, [SPECGRID] or [DIMENS]!'
        print('.....Done!')


        #Genetrate TOPS for cartesian grid if TOPS if not given
        if(self.grid_type=='Cartesian' and len(self.tops)==0):
            tops=np.zeros(self.n)
            for k in range(self.nz-1):
                for j in range(self.ny):
                    for i in range(self.nx):
                        ijk=cell_id(i,j,k,self.nx,self.ny)
                        ijk_next=cell_id(i,j,k+1,self.nx,self.ny)
                        tops[ijk_next] = tops[ijk] + self.dz[ijk]
            self.tops=tops


    def to_ecl(self, filename=None, keywords=None, one_file=False, return_string=False, save_file=True):
        
        
        list_str =[]
        key_added = []
        
        if self.grid_type == 'cartesian':

            if keywords is None:
                keywords = ['TOPS','DX','DY','DZ']
            else:
                assert isinstance(keywords,list)
            keywords_type = [i for i in keywords if i in ['TOPS','DX','DY','DZ']]
        else: 
            if keywords is None:
                keywords = ['COORD','ZCORN','SPECGRID']
            else:
                assert isinstance(keywords,list)
            keywords_type = [i for i in keywords if i in ['COORD','ZCORN']]
            
            if 'SPECGRID' in keywords:
                print('SPECGRID')
                list_str.append(f'SPECGRID\n {self.nx} {self.ny} {self.nz} 1 F /\n') 
                key_added.append('SPECGRID')
                
        if len(keywords_type)>0:
            for k in keywords_type:
                print(k)
                key_str = ""
                key_str += f'{k.upper()}\n'
                key_str += ' ' + ' '.join([str(v) + '\n' if (i+1)%10==0 else str(v) for i,v in enumerate(getattr(self,k.lower()))]) + '/\n'
                list_str.append(key_str)
                key_added.append(k)
                    
        
        keywords_spatial = [i for i in keywords if i not in ['SPECGRID','COORD','ZCORN','TOPS','DX','DY','DZ']]

        if all([bool(self.spatial_data),len(keywords_spatial)>0]):
            for key in keywords_spatial:
                print(key)
                key_str =""
                try:
                    key_str += key + '\n'
                    key_str += ' ' + ' '.join([str(v) + '\n' if (i+1)%10==0 else str(v) for i,v in enumerate(getattr(self,'spatial_data')[key])])  + '/\n'
                    list_str.append(key_str)   
                    key_added.append(key)                       
                except:
                    pass
        
        if save_file:            
            if one_file==True:
                if filename is None:
                    filename ='grid.GRDECL'
                try:
                    string = "".join(list_str)
                    with open(filename,'w') as text_file:
                        text_file.write(string)
                except Exception as e:
                    print(e)
                    pass
            
            else:
                if filename is None:
                    filename = '.'        
                filename = os.path.abspath(filename)
                print(filename)
                for i, key in enumerate(list_str):
                    
                    try:
                        with open(os.path.join(filename,key_added[i]+'.prop'),'w') as text_file:
                            text_file.write(key)
                    except Exception as e:
                        print(e)
                        pass
        
        if return_string:
            return "".join(list_str)
        
            
    
    def get_cell_id(self,i,j,k):
        """
        Get the cell Id given i,j,k indexes. 
            * ---  *  ---  *  --- *
            | 0,2  |  1,2  |  2,2 |  <- Cell id 6,7,8
            * ---  *  ---  *  --- *
            | 0,1  |  1,1  |  2,1 |  <- Cell id 3,4,5
            * ---  *  ---  *  --- *
            | 0,0  |  1,0  |  2,0 |  <- Cell id 0,1,2
            * ---  *  ---  * ---  *
        """
        c_id = cell_id(i,j,k,self.nx,self.ny)
        return c_id

    def get_cell_ijk(self,cell_id):
        """
        Get the cell indexes i,j,k given the cell id
            * ---  *  ---  *  --- *
            | 0,2  |  1,2  |  2,2 |  <- Cell id 6,7,8
            * ---  *  ---  *  --- *
            | 0,1  |  1,1  |  2,1 |  <- Cell id 3,4,5
            * ---  *  ---  *  --- *
            | 0,0  |  1,0  |  2,0 |  <- Cell id 0,1,2
            * ---  *  ---  * ---  *
        """
        i,j,k=cell_ijk(cell_id,self.nx,self.ny)
        return i,j,k

    def get_pillar(self,pillar_id:int):
        """
        Get the Top and Bottom coordinates of a pillar id
        """
        if self.grid_type == 'corner_point':
            id_top=[6*pillar_id+0,6*pillar_id+1,6*pillar_id+2]
            id_bottom=[6*pillar_id+3,6*pillar_id+4,6*pillar_id+5]
            top_point=np.array([self.coord[i] for i in id_top])
            bottom_point=np.array([self.coord[i] for i in id_bottom])
        else:
            raise ValueError('Pillar are only set in a Corner Point Grid')
        return np.array([top_point,bottom_point])

    def get_cell_pillars(self,i,j):
        """Obtain the four pillars (p0,p1,p2,p3) of a corner point cell
        The index of pillar
        
        3x3x1 system (2D X-Y plane)
        12--- 13  --- 14  ---15
        |      |       |      |  <- Cell 6,7,8
        8 ---  9  --- 10  ---11
        |      |       |      |  <- Cell 3,4,5
        4 ---  5  ---  6  --- 7
        |      |       |      |  <- Cell 0,1,2
        0 ---  1 ---   2 ---  3
        
        The pillars index for a grid follows below ordering (XY Plane)
        p2   p3
        *------*
        |      |
        |      |
        *------*
        p0   p1

        """
        if self.grid_type == 'corner_point':
            p0 = cell_id(i,j,0,self.nx+1,self.ny+1)
            p1 = cell_id(i+1,j,0,self.nx+1,self.ny+1)
            p2 = cell_id(i,j+1,0,self.nx+1,self.ny+1)
            p3 = cell_id(i+1,j+1,0,self.nx+1,self.ny+1)

            pls = [self.get_pillar(p0),self.get_pillar(p1),self.get_pillar(p2),self.get_pillar(p3)]
        else:
            raise ValueError('Pillar are only set in a Corner Point Grid')
        return np.array(pls)

    def get_vertices_id(self,i,j,k,order='GRD'):
        """
        Cartesian Grid

        Get the cell Id given i,j,k indexes. 
            13 --- 14  --- 15 --- 16
            | 0,2  |  1,2  |  2,2 |  <- Cell id 6,7,8
            9 ---  10  --- 11 --- 12
            | 0,1  |  1,1  |  2,1 |  <- Cell id 3,4,5
            5 ---  6  ---  7  --- 8
            | 0,0  |  1,0  |  2,0 |  <- Cell id 0,1,2
            1 ---  2  ---  3 ---  4

        Corner Point Grid

        3x3x1 system (2D X-Y plane)
        30---31,32---33,34---35
        |      |       |      |  <- Cell 6,7,8
        24---25,26---27,28---29
        18---19,20---21,22---23
        |      |       |      |  <- Cell 3,4,5
        12---13,14---15,16---17
        6 --- 7,8 --- 9,10---11
        |      |       |      |  <- Cell 0,1,2
        0 --- 1,2 --- 3,4 --- 5
        Node order convention for a 3D cell
         6----7
         -   -   <-Bottom Face
        4----5
          2----3
         -    -  <-Top Face
        0----1
        """ 
        if self.grid_type == 'cartesian':
            nx,ny=self.nx+1,self.ny+1
            p0=cell_id(i,j,k,nx,ny)
            p1=cell_id(i+1,j,k,nx,ny)
            p2=cell_id(i,j+1,k,nx,ny)
            p3=cell_id(i+1,j+1,k,nx,ny)

            p4=cell_id(i,j,k+1,nx,ny)
            p5=cell_id(i+1,j,k+1,nx,ny)
            p6=cell_id(i,j+1,k+1,nx,ny)
            p7=cell_id(i+1,j+1,k+1,nx,ny)

            if order == 'GRD':
                points = [p0,p1,p2,p3,p4,p5,p6,p7]
            elif order == 'VTK':
                points = [p4,p5,p7,p6,p0,p1,p3,p2]

            return np.array(points)

        if self.grid_type == 'corner_point':
            nx,ny=2*self.nx,2*self.ny
            p0=cell_id(2*i,2*j,2*k,nx,ny)
            p1=cell_id(2*i+1,2*j,2*k,nx,ny)
            p2=cell_id(2*i,2*j+1,2*k,nx,ny)
            p3=cell_id(2*i+1,2*j+1,2*k,nx,ny)

            p4=cell_id(2*i,2*j,2*k+1,nx,ny)
            p5=cell_id(2*i+1,2*j,2*k+1,nx,ny)
            p6=cell_id(2*i,2*j+1,2*k+1,nx,ny)
            p7=cell_id(2*i+1,2*j+1,2*k+1,nx,ny)

            if order == 'GRD':
                points = [p0,p1,p2,p3,p4,p5,p6,p7]
            elif order == 'VTK':
                points = [p4,p5,p7,p6,p0,p1,p3,p2]

            return np.array(points)


    def get_vertices_z(self,i,j,k):
        """
        Node order convention for a 3D cell
         6----7
         -   -   <-Bottom Face
        4----5
          2----3
         -    -  <-Top Face
        0----1
        """
        # Get the z coord for a cell
        if self.grid_type == 'corner_point':
            p = self.get_vertices_id(i,j,k)
            z = [self.zcorn[i] for i in p]
            return np.array(z)
        elif self.grid_type == 'cartesian':
            p= self.get_vertices_id(i,j,k)
            z = [self.cartesian_vertices_coord()[i,2] for i in p]
            return np.array(z)

        # Pending for cartessian grid

    def get_vertices_coords(self,i,j,k,order='GRD'):
        if self.grid_type == 'corner_point':
            coords=[]
            pillars = self.get_cell_pillars(i,j)
            cell_z =self.get_vertices_z(i,j,k)

            for i in range(8):
                p_id = i%4
                coords.append(interpolate_z_pillar(cell_z[i],pillars[p_id]))
            
            if order == 'GRD':
                v_coord = np.array(coords)
            elif order == 'VTK':
                v_coord = np.array(coords)[[4,5,7,6,0,1,3,2],:]

            return v_coord

        elif self.grid_type == 'cartesian':
            p= self.get_vertices_id(i,j,k)
            coords = [[self.cartesian_vertices_coord()[i,0],self.cartesian_vertices_coord()[i,1],self.cartesian_vertices_coord()[i,2]] for i in p]
            if order == 'GRD':
                v_coord = np.array(coords)
            elif order == 'VTK':
                v_coord = np.array(coords)[[4,5,7,6,0,1,3,2],:]

            return v_coord


    def get_vertices_face_z(self,i,j,k, face=None):
        """
         Get the Z coords for a cell
        
         6----7
         -   -   <-Bottom Face
        4----5
          2----3
         -    -  <-Top Face
        0----1   
        Follow getCornerPointCellIdx convention:
        X-, [0,2,4,6]
        X+, [1,3,5,7]
        Y-, [0,1,4,5]
        Y+, [2,3,6,7]
        Z+, [0,1,2,3]
        Z-, [4,5,6,7]
        """
        assert face is not None, 'A face must be choosen'
        points_id = self.get_vertices_id(i,j,k)
        if(face=="X-"): face_id=[points_id[0],points_id[2],points_id[4],points_id[6]]
        if(face=="X+"): face_id=[points_id[1],points_id[3],points_id[5],points_id[7]]
        if(face=="Y-"): face_id=[points_id[0],points_id[1],points_id[4],points_id[5]]
        if(face=="Y+"): face_id=[points_id[2],points_id[3],points_id[6],points_id[7]]
        if(face=="Z-"): face_id=[points_id[4],points_id[5],points_id[6],points_id[7]]
        if(face=="Z+"): face_id=[points_id[0],points_id[1],points_id[2],points_id[3]]

        if self.grid_type == 'cartesian':
            z_face = [self.cartesian_vertices_coord()[i,2] for i in face_id]
        elif self.grid_type == 'corner_point':
            z_face = [self.zcorn[i] for i in face_id]

        return np.array(z_face)

    def get_vertices_face_coords(self,i,j,k, face=None):
        """
         Get the Z coords for a cell
        
         6----7
         -   -   <-Bottom Face
        4----5
          2----3
         -    -  <-Top Face
        0----1   
        Follow getCornerPointCellIdx convention:
        X-, [0,2,4,6]
        X+, [1,3,5,7]
        Y-, [0,1,4,5]
        Y+, [2,3,6,7]
        Z+, [0,1,2,3]
        Z-, [4,5,6,7]
        """
        assert face is not None, 'A face must be choosen'
        points_id = self.get_vertices_id(i,j,k)
        if (face=="X-"): 
            face_id=[points_id[0],points_id[2],points_id[4],points_id[6]]
            ind = [0,2,4,6]
        elif (face=="X+"): 
            face_id=[points_id[1],points_id[3],points_id[5],points_id[7]]
            ind = [1,3,5,7]
        elif (face=="Y-"): 
            face_id=[points_id[0],points_id[1],points_id[4],points_id[5]]
            ind = [0,1,4,5]
        elif (face=="Y+"): 
            face_id=[points_id[2],points_id[3],points_id[6],points_id[7]]
            ind = [2,3,6,7]
        elif (face=="Z-"): 
            face_id=[points_id[4],points_id[5],points_id[6],points_id[7]]
            ind = [4,5,6,7]
        elif (face=="Z+"): 
            face_id=[points_id[0],points_id[1],points_id[2],points_id[3]]
            ind = [0,1,2,3]

        if self.grid_type == 'cartesian':
            z_face = [self.cartesian_vertices_coord()[i,:] for i in face_id]
        elif self.grid_type == 'corner_point':
            v_cord = self.get_vertices_coords(i,j,k)
            z_face = v_cord[ind,:]

        return np.array(z_face)


    def get_center_coord(self,i,j,k):
        """
        Get the cell Id given i,j,k indexes. 
            * ---  *  ---  *  --- *
            | 0,2  |  1,2  |  2,2 |  <- Cell id 6,7,8
            * ---  *  ---  *  --- *
            | 0,1  |  1,1  |  2,1 |  <- Cell id 3,4,5
            * ---  *  ---  *  --- *
            | 0,0  |  1,0  |  2,0 |  <- Cell id 0,1,2
            * ---  *  ---  * ---  *
        """
        cid = self.get_cell_id(i,j,k)
        if self.grid_type == 'cartesian':
            center = self.cartesian_center_point_coord()[cid,:]
        elif self.grid_type == 'corner_point':
            points = self.get_vertices_coords(i,j,k)
            center = points.mean(axis=0)
        
        return center

    def get_vtk(self):
        """
        Get the pyvista Object
        https://docs.pyvista.org/examples/00-load/create-unstructured-surface.html#sphx-glr-examples-00-load-create-unstructured-surface-py
        """
 
        #Identify the cell data connections
        offset = np.arange(0,9*self.n,step=9)

        points = np.zeros((self.n*8,3))
        #Cells
        for k in range(self.nz):
            for j in range(self.ny):
                for i in range(self.nx):
                    c_idx = self.get_cell_id(i,j,k)
                    #cells_array[c_idx,:] = self.get_vertices_id(i,j,k, order='VTK')

                    ind_from = 8*c_idx
                    ind_to = 8*(c_idx+1)
                    points[ind_from:ind_to,:] = self.get_vertices_coords(i,j,k, order = 'VTK')

        # Make a vector of shape self.n, make 2D and append to cell array then flatten C order
        cell_array = np.arange(self.n*8).reshape((self.n,8))
        cells = np.append(np.full(self.n,8).reshape((self.n,1)),cell_array,1).flatten()

        # cell type array. Contains the cell type of each cell
        cell_type = np.array([vtk.VTK_HEXAHEDRON]*self.n)

        grid = pv.UnstructuredGrid(offset, cells, cell_type, points)

        if self.spatial_data is not None:
            for i in self.spatial_data.items():
                grid.cell_arrays[i[0]] = i[1]

        return grid







        


                






    
    

    

    




