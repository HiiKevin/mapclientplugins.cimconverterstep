from datetime import time
import numpy as np
import csv
import os
import glob
import natsort
import pandas as pd
import re
import warnings
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, art3d


class Subdivision_Surface():
    '''
    This class create a surface from the control mesh, based on Catmull
    Clark subdivision surface method
        ...
        Attributes
        ----------
        etPos : array-like
            Matrix containing model after subdivision
        matrix : array-like
            Subdivision matrix (from coarse mesh to 3D model)  
        numNodes : int
            Number of nodes (coarse mesh)  
        numElements : int
            Number of elements 
        ETIndices : array-like
            Faces (subdivided model)
        control_mesh : array-like
            Coarse mesh 
        faceCM : array-like
            Coarse mesh connections
        faceCMsorted : array-like
            Faces (coarse mesh/template)
        etVertexStartEnd : array-like
            Start and End indices 
        CMStartEnd : array-like
            Start and End indices 
        SurfaceStartEnd : array-like
            Start and End indices
                    
        Methods
        -------
        Plot(SurfaceType,time_frame = 0)
            Plots the points of the model at a specified time frame on a 3d axis for visualisation.
        PlotCM(SurfaceType,time_frame = 0)
            Plots the points of the coarse mesh at a specified time frame on a 3d axis for visualisation.
    '''   
    etPos = None                # model after subdivision 
    matrix = None               # subdivision matrix (from coarse mesh to 3D model)        
    numNodes = 388              # number of nodes (coarse mesh)                 
    numElements = 187;          # number of elements 
    ETIndices = None            # Faces (subdivided model)
    control_mesh = None         # Coarse mesh 
    faceCM = None               # Coarse mesh connections 
    faceCMsorted = None         # Faces (coarse mesh/template)
    
    def __init__(self,type,FilePath):
        """
        Inputs
        ------
        type: str
            initialisation type
        FilePath: path-like
            directory where the model folder(s) are stored
        """

        dir_path=os.path.dirname(os.path.realpath(__file__))
        # load in matrices from text
        self.matrix = np.loadtxt(dir_path+'\\subdivision_matrix.txt') 
        self.ETIndices = np.loadtxt(dir_path+'\\ETIndices.txt')
        self.faceCM = np.loadtxt(dir_path+'\\ETIndices_control_mesh.txt')
        # initialise matrix indices
        self.CMStartEnd = np.array([[1,96],   # LV
                                    [97,211], # RV
                                    [212,354] # Epi
                                    ])-1      # convert indexing from 1 index to 0 index
        self.etVertexStartEnd = np.array([[1,1500],         # LV (vertices 1 to 1500 belong to the LV chamber...)
                                          [1501,2165],      # RVS
                                          [2166,3224],      # RVFW
                                          [3225,5582],      # Epi
                                          [5583,5631],      # Mitral valve
                                          [5632,5656],      # Aortic valve
                                          [5657,5697],      # Tricuspid valve
                                          [5698,5730],      # Pulmonary valve
                                          [5731,5810]       # RV insert
                                          ])-1              # convert indexing from 1 index to 0 index
        self.SurfaceStartEnd = np.array([[1,3072],          # LV (vertices 1 to 3072 belong to the LV chamber...)
                                         [3073,4480],       # RVS
                                         [4481,6752],       # RVFW
                                         [6753,11616],      # Epi
                                         [11617,11664],     # Mitral valve
                                         [11665,11688],     # Aortic valve
                                         [11689,11728],     # Tricuspid valve
                                         [11729,11760]      # Pulmonary valve
                                         ])-1               # convert indexing from 1 index to 0 index

        if type == 'InitFromControlMesh':
            self.control_mesh[:,:,1] = np.loadtxt(FilePath); 
            self.etPos[:,:,1] = np.matmul(self.matrix,self.control_mesh)
        elif type =='InitFromCIM':
            if  len(glob.glob(os.path.abspath(FilePath)+"\\*model\\*")) > 0: # ensure model folder(s) exists
                # extract .model files and sort numerically
                files = glob.glob(os.path.abspath(FilePath)+"\\*model\\*.model") 
                files = natsort.natsorted(files)
                # read each model file and store data
                for i in range(0,len(files)):
                    P = ReadFromCIMRVLVModel(os.path.join(os.path.abspath(FilePath),files[i]))
                    if self.control_mesh is not None:
                        self.control_mesh=np.concatenate((self.control_mesh, np.array([[P.XParameters,P.YParameters,P.ZParameters]]).T),axis=2)
                    else:
                        self.control_mesh=np.array([[P.XParameters,P.YParameters,P.ZParameters]]).T     # initialise matrix on first time
                    if self.etPos is not None:    
                        self.etPos=np.concatenate((self.etPos,np.matmul(self.matrix,self.control_mesh[:,:,i])[:,:,np.newaxis]),axis=2)
                    else:
                        self.etPos = np.matmul(self.matrix,self.control_mesh[:,:,i])[:,:,np.newaxis]    # initialise matrix on first time
            else:
                raise Exception('No CIM model folder found')
        else:
            raise Exception('please select your initialisation: InitFromCIM if you are dealing with CIM data or InitFromControlMesh if you only have a coarse mesh')
        
    def Plot(self,SurfaceType,time_frame = 0):
        """Plots the points of the model at a specified time frame on a 3d axis for visualisation.

        Inputs
        ------
        SurfaceType: str
            Type of surface that is to be plotted
        time_frame: int, optional
            Index of which time frame is to be plotted (default is 0)
        """

        if SurfaceType=='endo':
            faces_lv=self.ETIndices[range(self.SurfaceStartEnd[0,0],self.SurfaceStartEnd[0,1]+1),:]
            faces_rv=self.ETIndices[range(self.SurfaceStartEnd[1,0],self.SurfaceStartEnd[2,1]+1),:]
            x,y,z=self.etPos[:,:,time_frame].T
            ax = plt.figure().add_subplot(projection='3d')
            ax.set_box_aspect((np.max(x)-np.min(x), np.max(y)-np.min(y), np.max(z)-np.min(z)))
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            plt.title('Endocardium at t = '+str(time_frame))
            ax.plot_trisurf(x, y, z, triangles=faces_lv-1,color='green', alpha=0.75,edgecolors='black',linewidth=0.2)
            ax.plot_trisurf(x, y, z, triangles=faces_rv-1,color='blue', alpha=0.75,edgecolors='black',linewidth=0.2)
            plt.show()
        
        if SurfaceType=='LV':
            faces_lv=self.ETIndices[range(self.SurfaceStartEnd[0,0],self.SurfaceStartEnd[0,1]+1),:]
            x,y,z=self.etPos[:,:,time_frame].T
            ax = plt.figure().add_subplot(projection='3d')
            ax.set_box_aspect((np.max(x)-np.min(x), np.max(y)-np.min(y), np.max(z)-np.min(z)))
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            plt.title('LV at t = '+str(time_frame))
            ax.plot_trisurf(x, y, z, triangles=faces_lv-1,color='green', alpha=0.75,edgecolors='black',linewidth=0.2)
            plt.show()
        
        if SurfaceType=='RV':
            faces_rv=self.ETIndices[range(self.SurfaceStartEnd[1,0],self.SurfaceStartEnd[2,1]+1),:]
            x,y,z=self.etPos[:,:,time_frame].T
            ax = plt.figure().add_subplot(projection='3d')
            ax.set_box_aspect((np.max(x)-np.min(x), np.max(y)-np.min(y), np.max(z)-np.min(z)))
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            plt.title('RV at t = '+str(time_frame))
            ax.plot_trisurf(x, y, z, triangles=faces_rv-1,color='blue', alpha=0.75,edgecolors='black',linewidth=0.2)
            plt.show()

        if SurfaceType=='epi':
            faces_epi=self.ETIndices[range(self.SurfaceStartEnd[3,0],self.SurfaceStartEnd[3,1]+1),:]
            x,y,z=self.etPos[:,:,time_frame].T
            ax = plt.figure().add_subplot(projection='3d')
            ax.set_box_aspect((np.max(x)-np.min(x), np.max(y)-np.min(y), np.max(z)-np.min(z)))
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            plt.title('Epicardium at t = '+str(time_frame))
            ax.plot_trisurf(x, y, z, triangles=faces_epi-1,color='red', alpha=0.75,edgecolors='black',linewidth=0.2)
            plt.show()

        if SurfaceType=='all':
            faces_lv=self.ETIndices[range(self.SurfaceStartEnd[0,0],self.SurfaceStartEnd[0,1]+1),:]
            faces_rv=self.ETIndices[range(self.SurfaceStartEnd[1,0],self.SurfaceStartEnd[2,1]+1),:]
            faces_epi=self.ETIndices[range(self.SurfaceStartEnd[3,0],self.SurfaceStartEnd[3,1]+1),:]
            x,y,z=self.etPos[:,:,time_frame].T
            ax = plt.figure().add_subplot(projection='3d')
            ax.set_box_aspect((np.max(x)-np.min(x), np.max(y)-np.min(y), np.max(z)-np.min(z)))
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            plt.title('Biventricular model at t = '+str(time_frame))
            ax.plot_trisurf(x, y, z, triangles=faces_lv-1,color='green', alpha=0.9,edgecolors='black',linewidth=0.2)
            ax.plot_trisurf(x, y, z, triangles=faces_rv-1,color='blue', alpha=0.9,edgecolors='black',linewidth=0.2)
            ax.plot_trisurf(x, y, z, triangles=faces_epi-1,color='red', alpha=0.4,edgecolors='black',linewidth=0.2)
            plt.show()
    
    def PlotCM(self,SurfaceType,time_frame = 0):
        """Plots the points of the coarse mesh at a specified time frame on a 3d axis for visualisation.

        Inputs
        ------
        SurfaceType: str
            Type of surface that is to be plotted
        time_frame: int, optional
            Index of which time frame is to be plotted (default is 0)
        """

        if SurfaceType=='all':
            faces_lv=self.faceCM[range(self.CMStartEnd[0,0],self.CMStartEnd[0,1]+1),:]
            faces_rv=self.faceCM[range(self.CMStartEnd[1,0],self.CMStartEnd[1,1]+1),:]
            faces_epi=self.faceCM[range(self.CMStartEnd[2,0],self.CMStartEnd[2,1]+1),:]
            v=self.control_mesh[:,:,time_frame].T
            ax = plt.figure().add_subplot(projection='3d')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            ax.set_xlim3d(min(v[0,:])-1, max(v[0,:])+1)
            ax.set_ylim3d(min(v[1,:])-1, max(v[1,:])+1)
            ax.set_zlim3d(min(v[2,:])-1, max(v[2,:])+1)
            plt.title('Coarse Mesh at t = '+str(time_frame))
            pc1 = art3d.Poly3DCollection(v[:,(faces_lv-1).astype(int).T].T, facecolors='green', edgecolor="black",alpha=0.2)
            pc2 = art3d.Poly3DCollection(v[:,(faces_rv-1).astype(int).T].T, facecolors='blue', edgecolor="black",alpha=0.2)
            pc3 = art3d.Poly3DCollection(v[:,(faces_epi-1).astype(int).T].T, facecolors='red', edgecolor="black",alpha=0.2)
            ax.add_collection(pc1)
            ax.add_collection(pc2)
            ax.add_collection(pc3)
            plt.show()

        if SurfaceType=='LV':
            faces_lv=self.faceCM[range(self.CMStartEnd[0,0],self.CMStartEnd[0,1]+1),:]
            v=self.control_mesh[:,:,time_frame].T
            ax = plt.figure().add_subplot(projection='3d')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            ax.set_xlim3d(min(v[0,:])-1, max(v[0,:])+1)
            ax.set_ylim3d(min(v[1,:])-1, max(v[1,:])+1)
            ax.set_zlim3d(min(v[2,:])-1, max(v[2,:])+1)
            plt.title('Coarse LV at t = '+str(time_frame))
            pc = art3d.Poly3DCollection(v[:,(faces_lv-1).astype(int).T].T, facecolors='green', edgecolor="black",alpha=0.2)
            ax.add_collection(pc)
            plt.show()

        if SurfaceType=='RV':
            faces_rv=self.faceCM[range(self.CMStartEnd[1,0],self.CMStartEnd[1,1]+1),:]
            v=self.control_mesh[:,:,time_frame].T
            ax = plt.figure().add_subplot(projection='3d')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            ax.set_xlim3d(min(v[0,:])-1, max(v[0,:])+1)
            ax.set_ylim3d(min(v[1,:])-1, max(v[1,:])+1)
            ax.set_zlim3d(min(v[2,:])-1, max(v[2,:])+1)
            plt.title('Coarse RV at t = '+str(time_frame))
            pc = art3d.Poly3DCollection(v[:,(faces_rv-1).astype(int).T].T, facecolors='blue', edgecolor="black",alpha=0.2)
            ax.add_collection(pc)
            plt.show()

        if SurfaceType=='endo':
            faces_endo=self.faceCM[range(self.CMStartEnd[0,0],self.CMStartEnd[1,1]+1),:]
            v=self.control_mesh[:,:,time_frame].T
            ax = plt.figure().add_subplot(projection='3d')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            ax.set_xlim3d(min(v[0,:])-1, max(v[0,:])+1)
            ax.set_ylim3d(min(v[1,:])-1, max(v[1,:])+1)
            ax.set_zlim3d(min(v[2,:])-1, max(v[2,:])+1)
            plt.title('Coarse endocardium at t = '+str(time_frame))
            pc = art3d.Poly3DCollection(v[:,(faces_endo-1).astype(int).T].T, facecolors='red', edgecolor="black",alpha=0.2)
            ax.add_collection(pc)
            plt.show()

        if SurfaceType=='epi':
            faces_epi=self.faceCM[range(self.CMStartEnd[2,0],self.CMStartEnd[2,1]+1),:]
            v=self.control_mesh[:,:,time_frame].T
            ax = plt.figure().add_subplot(projection='3d')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            ax.set_xlim3d(min(v[0,:])-1, max(v[0,:])+1)
            ax.set_ylim3d(min(v[1,:])-1, max(v[1,:])+1)
            ax.set_zlim3d(min(v[2,:])-1, max(v[2,:])+1)
            plt.title('Coarse epicardium at t = '+str(time_frame))
            pc = art3d.Poly3DCollection(v[:,(faces_epi-1).astype(int).T].T, facecolors='red', edgecolor="black",alpha=0.2)
            ax.add_collection(pc)
            plt.show()


class ReadFromCIMRVLVModel():
    '''
    This class creates an RVLV model by reading from from a CIM RVLV file
        ...
        Attributes
        ----------
        ScaleFactor : float
            Scale factor
        ModelToMagnetTransform : array-like
            Model to magnet transform matrix
        XParameters : array-like
            X coordinates  
        YParameters : array-like
            Y coordinates  
        ZParameters : array-like
            Z coordinates  
    '''  

    def __init__(self,cim_file):
        """
        Inputs
        ------
        cim_file: path-like
            directory of CIM (.model) file to be read
        """

        # initialise dictionary for parameter types
        params = ['scaleFactor'.casefold(),
                  'ModelToMagnetTransform:ScaleFactor'.casefold(),
                  'ModelToMagnetTransform'.casefold(),'XParameters'.casefold(),
                  'YParameters'.casefold(),'ZParameters'.casefold()]
        df=pd.DataFrame({"type": pd.Series(['scalar','scalar','invIJK','vector','vector','vector'],index=params),
                         "value": pd.Series([[],[],[],[],[],[]], index=params)
                        })
        # ensure files can be read
        assert os.path.isfile(cim_file), 'Cannot read file ' + cim_file
        # read each .model file
        with open(cim_file, "r") as file:
            line = file.readline()
            while len(line) != 0:
                line=line.rstrip()
                # skip empty line
                if len(line) == 0:
                    line = file.readline()
                    continue
                # ignore this one
                if line == 'ModelToMagnetTransform:':
                    line = file.readline()
                    continue
                # remove last
                if line[-1] == ':':
                    line=line[0:-1]
                # find which parameter it is
                assert line.casefold() in df.index, 'unhandled parameter type '+line
                param=line.casefold()
                if df['type'][param] == 'scalar':
                    df['value'][param].append(float(file.readline().rstrip()))
                elif df['type'][param] == 'vector':
                    nelmt=int(file.readline().rstrip())
                    for i in range(0,nelmt):
                        df['value'][param].append(float(file.readline().rstrip()))
                elif df['type'][param] == 'invIJK':
                    for i in range(0,4):
                        data=[float(i) for i in re.findall('([0-9e.\-]+)i\s+([0-9e.\-]+)j\s+([0-9e.\-]+)k',file.readline().rstrip())[0]]
                        df['value'][param].append(data)
                else:
                    warnings.warn('Unknown parameter type of.' + param + 'Skip.\n',)
                    continue
                line=file.readline()
        file.close() 
        # store data in appropriate arrays
        self.scaleFactor=df['value']['scaleFactor'.casefold()]
        self.ModelToMagnetTransform=np.array(df['value']['ModelToMagnetTransform'.casefold()]).T
        self.XParameters=np.array(df['value']['XParameters'.casefold()])
        self.YParameters=np.array(df['value']['YParameters'.casefold()])
        self.ZParameters=np.array(df['value']['ZParameters'.casefold()])
        

def main():
    """This script reads the cardiohance_082 .model files to create a Subdivision_Surface
    model and writes the xyz coordinates of the model to csv files in the OutputPath.
    """

    FilePath = '.\\input\\cardiohance_082'
    OutputPath = '.\\output_python'
    model = Subdivision_Surface('InitFromCIM',FilePath)
    time_frame = np.shape(model.etPos)[2]
    # write to csv (disabled for debugging)
    if False:
        for i in range(time_frame):
                LVendo = model.etPos[range(model.etVertexStartEnd[0,0],model.etVertexStartEnd[0,1]+1),:,i]
                RV_S = model.etPos[range(model.etVertexStartEnd[1,0],model.etVertexStartEnd[1,1]+1),:,i]
                RV_FW = model.etPos[range(model.etVertexStartEnd[2,0],model.etVertexStartEnd[2,1]+1),:,i]
                Epi = model.etPos[range(model.etVertexStartEnd[3,0],model.etVertexStartEnd[3,1]+1),:,i]
                MV = model.etPos[range(model.etVertexStartEnd[4,0],model.etVertexStartEnd[4,1]+1),:,i]
                AV = model.etPos[range(model.etVertexStartEnd[5,0],model.etVertexStartEnd[5,1]+1),:,i]
                TV = model.etPos[range(model.etVertexStartEnd[6,0],model.etVertexStartEnd[6,1]+1),:,i]
                PV = model.etPos[range(model.etVertexStartEnd[7,0],model.etVertexStartEnd[7,1]+1),:,i]
                np.savetxt(os.path.abspath(OutputPath)+'\\'+'LVendo_'+str(i+1)+'.csv',LVendo,delimiter=',')
                np.savetxt(os.path.abspath(OutputPath)+'\\'+'RV_septum_'+str(i+1)+'.csv',RV_S,delimiter=',')
                np.savetxt(os.path.abspath(OutputPath)+'\\'+'RV_freewall_'+str(i+1)+'.csv',RV_FW,delimiter=',')
                np.savetxt(os.path.abspath(OutputPath)+'\\'+'Epi_'+str(i+1)+'.csv',Epi,delimiter=',')
                np.savetxt(os.path.abspath(OutputPath)+'\\'+'MV_'+str(i+1)+'.csv',MV[:-1,:],delimiter=',')
                np.savetxt(os.path.abspath(OutputPath)+'\\'+'AV_'+str(i+1)+'.csv',AV[:-1,:],delimiter=',')
                np.savetxt(os.path.abspath(OutputPath)+'\\'+'TV_'+str(i+1)+'.csv',TV[:-1,:],delimiter=',')
                np.savetxt(os.path.abspath(OutputPath)+'\\'+'PV_'+str(i+1)+'.csv',PV[:-1,:],delimiter=',')

if __name__ == "__main__":
    main()