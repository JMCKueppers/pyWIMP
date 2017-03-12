from ROOT import *
from pyWIMP_Halo.utilities.hist_tools import *

def maxwell_boltzmann(v_sub_E=235, v_sub_0=220, v_sub_esc=544):
  if v_sub_esc == 0:
    maxwell_boltzmann = TF1("Maxwell-Boltzmann", "TMath::Power(x,2)*TMath::Exp(-TMath::Power((x+[0]),2)/TMath::Power([1],2))", 0, 1000)
    maxwell_boltzmann.SetParameter(0, v_sub_E)
    maxwell_boltzmann.SetParameter(1, v_sub_0)
  else:
    maxwell_boltzmann = TF1("Maxwell-Boltzmann", "(x<=[2])*TMath::Power(x,2)*TMath::Exp(-TMath::Power((x+[0]),2)/TMath::Power([1],2))", 0, 1000)
    maxwell_boltzmann.SetParameter(0, v_sub_E)
    maxwell_boltzmann.SetParameter(1, v_sub_0)
    maxwell_boltzmann.SetParameter(2, v_sub_esc)
  
  export_tf1(maxwell_boltzmann, "%s/../../input_data/MaxwellBoltzmann.dat" % os.path.dirname(os.path.realpath(__file__)), 10000)
 
 
def rectangle(v_sub_start=0, v_sub_esc=544):
  if v_sub_esc == 0:
    rectangle = TF1("Rectangle", "(x>=[0])", 0, 1000)
    rectangle.SetParameter(0, v_sub_start)
  else:
    rectangle = TF1("Rectangle", "(x>=[0])*(x<=[1])", 0, 1000)
    rectangle.SetParameter(0, v_sub_start)
    rectangle.SetParameter(1, v_sub_esc)
  
  export_tf1(rectangle, "%s/../../input_data/Rectangle.dat" % os.path.dirname(os.path.realpath(__file__)), 10000)
 

def triangle(v_sub_start=0, v_sub_center=220, v_sub_stop=544, v_sub_esc=544):
  if v_sub_esc == 0:
    triangle = TF1("Triangle", "(x>=[0])*(x<=[1])*(1./([1]-[0]))*(x-[0])-(x>[1])*(x<=[2])*(1./([2]-[1]))*(x-[2])", 0, 1000)
    triangle.SetParameter(0, v_sub_start)
    triangle.SetParameter(1, v_sub_center)
    triangle.SetParameter(2, v_sub_stop)
  else:
    triangle = TF1("Triangle", "(x<=[3])*((x>=[0])*(x<=[1])*(1./([1]-[0]))*(x-[0])-(x>[1])*(x<=[2])*(1./([2]-[1]))*(x-[2]))", 0, 1000)
    triangle.SetParameter(0, v_sub_start)
    triangle.SetParameter(1, v_sub_center)
    triangle.SetParameter(2, v_sub_stop)
    triangle.SetParameter(3, v_sub_esc)

  export_tf1(triangle, "%s/../../input_data/Triangle.dat" % os.path.dirname(os.path.realpath(__file__)), 10000)


def circle(v_sub_center=220, radius=220, v_sub_esc=544):
  if v_sub_esc == 0:
    circle = TF1("Circle", "(x>=0)*TMath::Sqrt(TMath::Power([1],2)-TMath::Power(x-[0],2))", 0, 1000)
    circle.SetParameter(0, v_sub_center)
    circle.SetParameter(1, radius)
  else:
    circle = TF1("Circle", "(x>=0)*(x<=[2])*TMath::Sqrt(TMath::Power([1],2)-TMath::Power(x-[0],2))", 0, 1000)
    circle.SetParameter(0, v_sub_center)
    circle.SetParameter(1, radius)
    circle.SetParameter(2, v_sub_esc)
  
  export_tf1(circle, "%s/../../input_data/Circle.dat" % os.path.dirname(os.path.realpath(__file__)), 10000)


def gauss(v_sub_center=220, width=(544-220)/3, v_sub_esc=544):
  if v_sub_esc == 0:
    gauss = TF1("Gauss", "TMath::Exp(-0.5*TMath::Power((x-[0])/[1],2))", 0, 1000)
    gauss.SetParameter(0, v_sub_center)
    gauss.SetParameter(1, width)
  else:
    gauss = TF1("Gauss", "(x<=[2])*TMath::Exp(-0.5*TMath::Power((x-[0])/[1],2))", 0, 1000)
    gauss.SetParameter(0, v_sub_center)
    gauss.SetParameter(1, width)
    gauss.SetParameter(2, v_sub_esc)
  
  export_tf1(gauss, "%s/../../input_data/Gauss.dat" % os.path.dirname(os.path.realpath(__file__)), 10000)


def parabola(v_sub_center=220, width=(544-220)/3, v_sub_esc=544):
  if v_sub_esc == 0:
    parabola = TF1("Parabola", "TMath::Power([1],2)-TMath::Power(x-[0],2)", 0, 1000)
    parabola.SetParameter(0, v_sub_center)
    parabola.SetParameter(1, width)
  else:
    parabola = TF1("Parabola", "(x<=[2])*(TMath::Power([1],2)-TMath::Power(x-[0],2))", 0, 1000)
    parabola.SetParameter(0, v_sub_center)
    parabola.SetParameter(1, width)
    parabola.SetParameter(2, v_sub_esc)
  
  export_tf1(parabola, "%s/../../input_data/Parabola.dat" % os.path.dirname(os.path.realpath(__file__)), 10000)
