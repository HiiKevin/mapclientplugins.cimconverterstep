
"""
MAP Client Plugin Step
"""
import json
import os
import numpy as np
from mapclientplugins.cimconverterstep.SurfaceExtractor import Subdivision_Surface
from PySide2 import QtGui

from mapclient.mountpoints.workflowstep import WorkflowStepMountPoint
from mapclientplugins.cimconverterstep.configuredialog import ConfigureDialog


class CIMConverterStep(WorkflowStepMountPoint):
    """
    Skeleton step which is intended to be a helpful starting point
    for new steps.
    """

    def __init__(self, location):
        super(CIMConverterStep, self).__init__('CIM Converter', location)
        self._configured = False # A step cannot be executed until it has been configured.
        self._category = 'Source'
        # Add any other initialisation code here:
        self._icon =  QtGui.QImage(':/cimconverterstep/images/data-source.png')
        # Ports:
        self.addPort(('http://physiomeproject.org/workflow/1.0/rdf-schema#port',
                      'http://physiomeproject.org/workflow/1.0/rdf-schema#provides',
                      'http://physiomeproject.org/workflow/1.0/rdf-schema#file_location'))
        # Port data:
        self._portData0 = None # http://physiomeproject.org/workflow/1.0/rdf-schema#file_location
        # Config:
        self._config = {}
        self._config['identifier'] = ''

    def execute(self):
        """
        Add your code here that will kick off the execution of the step.
        Make sure you call the _doneExecution() method when finished.  This method
        may be connected up to a button in a widget for example.
        """
        # Put your execute step code here before calling the '_doneExecution' method.
        input_path = os.path.abspath(os.path.join(self._location, self._config['file']))
        output_path=os.path.join(input_path+'\\csv')
        if not os.path.isdir(output_path):
            os.makedirs(output_path)
        model = Subdivision_Surface('InitFromCIM',input_path)
        time_frame = np.shape(model.etPos)[2]
        for i in range(time_frame):
            LVendo = model.etPos[range(model.etVertexStartEnd[0,0],model.etVertexStartEnd[0,1]+1),:,i]
            RV_S = model.etPos[range(model.etVertexStartEnd[1,0],model.etVertexStartEnd[1,1]+1),:,i]
            RV_FW = model.etPos[range(model.etVertexStartEnd[2,0],model.etVertexStartEnd[2,1]+1),:,i]
            Epi = model.etPos[range(model.etVertexStartEnd[3,0],model.etVertexStartEnd[3,1]+1),:,i]
            MV = model.etPos[range(model.etVertexStartEnd[4,0],model.etVertexStartEnd[4,1]+1),:,i]
            AV = model.etPos[range(model.etVertexStartEnd[5,0],model.etVertexStartEnd[5,1]+1),:,i]
            TV = model.etPos[range(model.etVertexStartEnd[6,0],model.etVertexStartEnd[6,1]+1),:,i]
            PV = model.etPos[range(model.etVertexStartEnd[7,0],model.etVertexStartEnd[7,1]+1),:,i]
            np.savetxt(os.path.abspath(output_path)+'\\'+'LVendo_'+str(i+1)+'.csv',LVendo,delimiter=',')
            np.savetxt(os.path.abspath(output_path)+'\\'+'RV_septum_'+str(i+1)+'.csv',RV_S,delimiter=',')
            np.savetxt(os.path.abspath(output_path)+'\\'+'RV_freewall_'+str(i+1)+'.csv',RV_FW,delimiter=',')
            np.savetxt(os.path.abspath(output_path)+'\\'+'Epi_'+str(i+1)+'.csv',Epi,delimiter=',')
            np.savetxt(os.path.abspath(output_path)+'\\'+'MV_'+str(i+1)+'.csv',MV[:-1,:],delimiter=',')
            np.savetxt(os.path.abspath(output_path)+'\\'+'AV_'+str(i+1)+'.csv',AV[:-1,:],delimiter=',')
            np.savetxt(os.path.abspath(output_path)+'\\'+'TV_'+str(i+1)+'.csv',TV[:-1,:],delimiter=',')
            np.savetxt(os.path.abspath(output_path)+'\\'+'PV_'+str(i+1)+'.csv',PV[:-1,:],delimiter=',')        
        self._portData0=output_path
        self._doneExecution()

    def getPortData(self, index):
        """
        Add your code here that will return the appropriate objects for this step.
        The index is the index of the port in the port list.  If there is only one
        provides port for this step then the index can be ignored.

        :param index: Index of the port to return.
        """
        return self._portData0 # http://physiomeproject.org/workflow/1.0/rdf-schema#file_location

    def configure(self):
        """
        This function will be called when the configure icon on the step is
        clicked.  It is appropriate to display a configuration dialog at this
        time.  If the conditions for the configuration of this step are complete
        then set:
            self._configured = True
        """
        dlg = ConfigureDialog(self._main_window)
        dlg.setWorkflowLocation(self._location)
        dlg.identifierOccursCount = self._identifierOccursCount
        dlg.setConfig(self._config)
        dlg.validate()
        dlg.setModal(True)

        if dlg.exec_():
            self._config = dlg.getConfig()

        self._configured = dlg.validate()
        self._configuredObserver()

    def getIdentifier(self):
        """
        The identifier is a string that must be unique within a workflow.
        """
        return self._config['identifier']

    def setIdentifier(self, identifier):
        """
        The framework will set the identifier for this step when it is loaded.
        """
        self._config['identifier'] = identifier

    def serialize(self):
        """
        Add code to serialize this step to string.  This method should
        implement the opposite of 'deserialize'.
        """
        return json.dumps(self._config, default=lambda o: o.__dict__, sort_keys=True, indent=4)

    def deserialize(self, string):
        """
        Add code to deserialize this step from string.  This method should
        implement the opposite of 'serialize'.

        :param string: JSON representation of the configuration in a string.
        """
        self._config.update(json.loads(string))

        d = ConfigureDialog()
        d.setWorkflowLocation(self._location)
        d.identifierOccursCount = self._identifierOccursCount
        d.setConfig(self._config)
        self._configured = d.validate()


