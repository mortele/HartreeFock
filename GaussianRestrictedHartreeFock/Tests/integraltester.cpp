#include "integraltester.h"
#include "Orbitals/gaussianprimitive.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using std::fabs;
using arma::vec;
using arma::mat;
using arma::zeros;
using std::cout;
using std::endl;
using std::setprecision;

void IntegralTester::setupExactOverlapMatrix() {
    m_exactOverlap(0,0) = 1.96870124321530;
    m_exactOverlap(0,1) = -0.0341172323898471;
    m_exactOverlap(0,2) = -0.804031477774766;
    m_exactOverlap(0,3) = 0.527226809908349;
    m_exactOverlap(0,4) = 2.95367374025867e-6;
    m_exactOverlap(0,5) = 0.000675323284399051;
    m_exactOverlap(0,6) = 0.0202562111841077;
    m_exactOverlap(0,7) = 0.000967087011747088;
    m_exactOverlap(0,8) = 2.13813050004227e-5;
    m_exactOverlap(0,9) = 1.40984128295943e-19;
    m_exactOverlap(1,0) = -0.0341172323898471;
    m_exactOverlap(1,1) = 0.387827062976667;
    m_exactOverlap(1,2) = -0.00176248709791666;
    m_exactOverlap(1,3) = 0.0122388517922682;
    m_exactOverlap(1,4) = 4.99655646048901e-11;
    m_exactOverlap(1,5) = 0.00734798949097483;
    m_exactOverlap(1,6) = 0.00138996712571738;
    m_exactOverlap(1,7) = -0.00127323177225232;
    m_exactOverlap(1,8) = 0.00101665428542184;
    m_exactOverlap(1,9) = -1.81927486326609e-5;
    m_exactOverlap(2,0) = -0.804031477774766;
    m_exactOverlap(2,1) = -0.00176248709791666;
    m_exactOverlap(2,2) = 9.98427851756932;
    m_exactOverlap(2,3) = -0.0606071391980689;
    m_exactOverlap(2,4) = 0.0185755368982808;
    m_exactOverlap(2,5) = 0.000961861089058550;
    m_exactOverlap(2,6) = 0.0383258102957808;
    m_exactOverlap(2,7) = -0.000135744841834612;
    m_exactOverlap(2,8) = -8.30819528805685e-5;
    m_exactOverlap(2,9) = 0.000882925086794951;
    m_exactOverlap(3,0) = 0.527226809908349;
    m_exactOverlap(3,1) = 0.0122388517922682;
    m_exactOverlap(3,2) = -0.0606071391980689;
    m_exactOverlap(3,3) = 0.640491765531426;
    m_exactOverlap(3,4) = -0.00289162677173892;
    m_exactOverlap(3,5) = 4.59938683989704e-5;
    m_exactOverlap(3,6) = 0.0224417428337129;
    m_exactOverlap(3,7) = 0.000155680871167730;
    m_exactOverlap(3,8) = 1.28786125356206e-6;
    m_exactOverlap(3,9) = 0.000399808605221194;
    m_exactOverlap(4,0) = 2.95367374025867e-6;
    m_exactOverlap(4,1) = 4.99655646048901e-11;
    m_exactOverlap(4,2) = 0.0185755368982808;
    m_exactOverlap(4,3) = -0.00289162677173892;
    m_exactOverlap(4,4) = 0.0233723135930304;
    m_exactOverlap(4,5) = 1.51663773622826e-24;
    m_exactOverlap(4,6) = 3.63749494189881e-12;
    m_exactOverlap(4,7) = 4.69426663197159e-20;
    m_exactOverlap(4,8) = 1.03888645793693e-24;
    m_exactOverlap(4,9) = 2.17123310831415e-13;
    m_exactOverlap(5,0) = 0.000675323284399051;
    m_exactOverlap(5,1) = 0.00734798949097483;
    m_exactOverlap(5,2) = 0.000961861089058550;
    m_exactOverlap(5,3) = 4.59938683989704e-5;
    m_exactOverlap(5,4) = 1.51663773622826e-24;
    m_exactOverlap(5,5) = 0.0172362086060265;
    m_exactOverlap(5,6) = 8.26196083281730e-7;
    m_exactOverlap(5,7) = -2.20707377259169e-5;
    m_exactOverlap(5,8) = 0.00100873345509930;
    m_exactOverlap(5,9) = -1.27038048996005e-7;
    m_exactOverlap(6,0) = 0.0202562111841077;
    m_exactOverlap(6,1) = 0.00138996712571738;
    m_exactOverlap(6,2) = 0.0383258102957808;
    m_exactOverlap(6,3) = 0.0224417428337129;
    m_exactOverlap(6,4) = 3.63749494189881e-12;
    m_exactOverlap(6,5) = 8.26196083281730e-7;
    m_exactOverlap(6,6) = 0.00703744735730407;
    m_exactOverlap(6,7) = -1.54331741523399e-5;
    m_exactOverlap(6,8) = 2.01930457151725e-8;
    m_exactOverlap(6,9) = 5.69065367227716e-5;
    m_exactOverlap(7,0) = 0.000967087011747088;
    m_exactOverlap(7,1) = -0.00127323177225232;
    m_exactOverlap(7,2) = -0.000135744841834612;
    m_exactOverlap(7,3) = 0.000155680871167730;
    m_exactOverlap(7,4) = 4.69426663197159e-20;
    m_exactOverlap(7,5) = -2.20707377259169e-5;
    m_exactOverlap(7,6) = -1.54331741523399e-5;
    m_exactOverlap(7,7) = 0.00126285590311933;
    m_exactOverlap(7,8) = -7.49167731362069e-7;
    m_exactOverlap(7,9) = 5.91524141572004e-5;
    m_exactOverlap(8,0) = 2.13813050004227e-5;
    m_exactOverlap(8,1) = 0.00101665428542184;
    m_exactOverlap(8,2) = -8.30819528805685e-5;
    m_exactOverlap(8,3) = 1.28786125356206e-6;
    m_exactOverlap(8,4) = 1.03888645793693e-24;
    m_exactOverlap(8,5) = 0.00100873345509930;
    m_exactOverlap(8,6) = 2.01930457151725e-8;
    m_exactOverlap(8,7) = -7.49167731362069e-7;
    m_exactOverlap(8,8) = 0.0491199707835051;
    m_exactOverlap(8,9) = 2.13719828184536e-10;
    m_exactOverlap(9,0) = 1.40984128295943e-19;
    m_exactOverlap(9,1) = -1.81927486326609e-5;
    m_exactOverlap(9,2) = 0.000882925086794951;
    m_exactOverlap(9,3) = 0.000399808605221194;
    m_exactOverlap(9,4) = 2.17123310831415e-13;
    m_exactOverlap(9,5) = -1.27038048996005e-7;
    m_exactOverlap(9,6) = 5.69065367227716e-5;
    m_exactOverlap(9,7) = 5.91524141572004e-5;
    m_exactOverlap(9,8) = 2.13719828184536e-10;
    m_exactOverlap(9,9) = 0.000746311679740889;
}

void IntegralTester::setupExactKineticMatrix() {
    m_exactKinetic(0,0) = 2.95305186482296;
    m_exactKinetic(0,1) = -0.0558424051783208;
    m_exactKinetic(0,2) = -0.710211946758681;
    m_exactKinetic(0,3) = 0.825190282920817;
    m_exactKinetic(0,4) = -4.09145004078157e-5;
    m_exactKinetic(0,5) = -0.00354106288941588;
    m_exactKinetic(0,6) = -0.00101696409225696;
    m_exactKinetic(0,7) = 0.000318784100712331;
    m_exactKinetic(0,8) = -0.000193277296800795;
    m_exactKinetic(0,9) = -1.03461059186122e-16;
    m_exactKinetic(1,0) = -0.0558424051783208;
    m_exactKinetic(1,1) = 1.06652442318584;
    m_exactKinetic(1,2) = -0.00136283235772673;
    m_exactKinetic(1,3) = 0.00666163150823509;
    m_exactKinetic(1,4) = -1.38735470820744e-9;
    m_exactKinetic(1,5) = 0.0101718786886854;
    m_exactKinetic(1,6) = -0.00332120851907597;
    m_exactKinetic(1,7) = -0.00085909102921053;
    m_exactKinetic(1,8) = -0.00197612918039297;
    m_exactKinetic(1,9) = -0.000269901350336851;
    m_exactKinetic(2,0) = -0.710211946758697;
    m_exactKinetic(2,1) = -0.00136283235772661;
    m_exactKinetic(2,2) = 7.48820888817698;
    m_exactKinetic(2,3) = -0.0941675849397594;
    m_exactKinetic(2,4) = 0.00992457461917693;
    m_exactKinetic(2,5) = -0.00129538182364648;
    m_exactKinetic(2,6) = 0.0314870023540979;
    m_exactKinetic(2,7) = 8.81768102492484e-6;
    m_exactKinetic(2,8) = 0.000281024428562610;
    m_exactKinetic(2,9) = 0.00157572335530054;
    m_exactKinetic(3,0) = 0.825190282920816;
    m_exactKinetic(3,1) = 0.00666163150823466;
    m_exactKinetic(3,2) = -0.0941675849397592;
    m_exactKinetic(3,3) = 1.44110647244571;
    m_exactKinetic(3,4) = 0.00871686596713676;
    m_exactKinetic(3,5) = -0.000409048293848453;
    m_exactKinetic(3,6) = 0.0304280527464197;
    m_exactKinetic(3,7) = 2.05945927793234e-5;
    m_exactKinetic(3,8) = -1.43955534859779e-5;
    m_exactKinetic(3,9) = 0.00208775559366966;
    m_exactKinetic(4,0) = -4.09145004078137e-5;
    m_exactKinetic(4,1) = -1.38735470820744e-9;
    m_exactKinetic(4,2) = 0.00992457461936637;
    m_exactKinetic(4,3) = 0.00871686596713263;
    m_exactKinetic(4,4) = 0.111408028141933;
    m_exactKinetic(4,5) = -1.80773041407557e-22;
    m_exactKinetic(4,6) = -2.08223149360596e-10;
    m_exactKinetic(4,7) = -4.77727574575186e-18;
    m_exactKinetic(4,8) = -8.97002227888797e-23;
    m_exactKinetic(4,9) = -1.46810972477601e-11;
    m_exactKinetic(5,0) = -0.00354106288941604;
    m_exactKinetic(5,1) = 0.0101718786886830;
    m_exactKinetic(5,2) = -0.00129538182364563;
    m_exactKinetic(5,3) = -0.000409048293848444;
    m_exactKinetic(5,4) = -1.80773041407557e-22;
    m_exactKinetic(5,5) = 0.0896282847513165;
    m_exactKinetic(5,6) = -1.92843648148517e-5;
    m_exactKinetic(5,7) = -0.000172941244267707;
    m_exactKinetic(5,8) = -0.00462092936135390;
    m_exactKinetic(5,9) = 2.73212268997395e-6;
    m_exactKinetic(6,0) = -0.00101696409225621;
    m_exactKinetic(6,1) = -0.00332120851907591;
    m_exactKinetic(6,2) = 0.0314870023540918;
    m_exactKinetic(6,3) = 0.0304280527464185;
    m_exactKinetic(6,4) = -2.08223149360596e-10;
    m_exactKinetic(6,5) = -1.92843648148518e-5;
    m_exactKinetic(6,6) = 0.0472681880832255;
    m_exactKinetic(6,7) = -0.000206874924595757;
    m_exactKinetic(6,8) = -4.32226463603479e-7;
    m_exactKinetic(6,9) = 0.00105258842223646;
    m_exactKinetic(7,0) = 0.000318784100712428;
    m_exactKinetic(7,1) = -0.000859091029210779;
    m_exactKinetic(7,2) = 8.81768102490486e-6;
    m_exactKinetic(7,3) = 2.05945927793469e-5;
    m_exactKinetic(7,4) = -4.77727574575185e-18;
    m_exactKinetic(7,5) = -0.000172941244267711;
    m_exactKinetic(7,6) = -0.000206874924595625;
    m_exactKinetic(7,7) = 0.0163539839453953;
    m_exactKinetic(7,8) = 5.36473707571811e-6;
    m_exactKinetic(7,9) = 0.000246865723384732;
    m_exactKinetic(8,0) = -0.000193277296800797;
    m_exactKinetic(8,1) = -0.00197612918039266;
    m_exactKinetic(8,2) = 0.000281024428562512;
    m_exactKinetic(8,3) = -1.43955534859778e-5;
    m_exactKinetic(8,4) = -8.97002227888790e-23;
    m_exactKinetic(8,5) = -0.00462092936135679;
    m_exactKinetic(8,6) = -4.32226463603479e-7;
    m_exactKinetic(8,7) = 5.36473707571820e-6;
    m_exactKinetic(8,8) = 0.223495867064918;
    m_exactKinetic(8,9) = -4.82109324561387e-9;
    m_exactKinetic(9,0) = -1.50727578348677e-17;
    m_exactKinetic(9,1) = -0.000269901350336876;
    m_exactKinetic(9,2) = 0.00157572335530140;
    m_exactKinetic(9,3) = 0.00208775559366982;
    m_exactKinetic(9,4) = -1.46810972477601e-11;
    m_exactKinetic(9,5) = 2.73212268997395e-6;
    m_exactKinetic(9,6) = 0.00105258842223629;
    m_exactKinetic(9,7) = 0.000246865723384709;
    m_exactKinetic(9,8) = -4.82109324561561e-9;
    m_exactKinetic(9,9) = 0.0112319907800976;
}

void IntegralTester::setupExactElectronNucleusMatrix() {
    m_exactElectronNucleus(0,0) = 3.141592653396305;
    m_exactElectronNucleus(0,1) = 0.006074920002336657;
    m_exactElectronNucleus(0,2) = -0.23052941394817655;
    m_exactElectronNucleus(0,3) = 0.46132628910842943;
    m_exactElectronNucleus(0,4) = 7.268945590080572e-7;
    m_exactElectronNucleus(0,5) = 1.8434518400328907e-6;
    m_exactElectronNucleus(0,6) = 0.027310496584566385;
    m_exactElectronNucleus(0,7) = 0.0001854033428543244;
    m_exactElectronNucleus(0,8) = 5.778111260020206e-6;
    m_exactElectronNucleus(0,9) = -3.5809106251696886e-19;
    m_exactElectronNucleus(1,0) = -0.013094721102328675;
    m_exactElectronNucleus(1,1) = 0.43272621892355256;
    m_exactElectronNucleus(1,2) = 0.0030270630807894675;
    m_exactElectronNucleus(1,3) = 0.003442049393588565;
    m_exactElectronNucleus(1,4) = 1.4228902001479771e-11;
    m_exactElectronNucleus(1,5) = 0.003430101670771089;
    m_exactElectronNucleus(1,6) = 0.0006314855924455602;
    m_exactElectronNucleus(1,7) = -0.0001807195948772219;
    m_exactElectronNucleus(1,8) = 0.0002691063236345793;
    m_exactElectronNucleus(1,9) = -1.0462479019948222e-5;
    m_exactElectronNucleus(2,0) = -0.4510728322381398;
    m_exactElectronNucleus(2,1) = 0.010286345926581059;
    m_exactElectronNucleus(2,2) = 5.8177641728586185;
    m_exactElectronNucleus(2,3) = 0.04061181557716256;
    m_exactElectronNucleus(2,4) = 0.007531157186406296;
    m_exactElectronNucleus(2,5) = 0.00042066329652329707;
    m_exactElectronNucleus(2,6) = 0.008868246560207824;
    m_exactElectronNucleus(2,7) = 2.1136695773265444e-5;
    m_exactElectronNucleus(2,8) = 7.92402055095499e-5;
    m_exactElectronNucleus(2,9) = 0.0002839888511347502;
    m_exactElectronNucleus(3,0) = 0.4880275464523432;
    m_exactElectronNucleus(3,1) = 0.00264818690665339;
    m_exactElectronNucleus(3,2) = -0.07678460604402675;
    m_exactElectronNucleus(3,3) = 0.5166079677373133;
    m_exactElectronNucleus(3,4) = -0.0006552403098904724;
    m_exactElectronNucleus(3,5) = 3.9867724278037694e-5;
    m_exactElectronNucleus(3,6) = 0.017774616398394436;
    m_exactElectronNucleus(3,7) = 0.0003086285665016171;
    m_exactElectronNucleus(3,8) = 5.559844717322346e-7;
    m_exactElectronNucleus(3,9) = 0.0003565006556826383;
    m_exactElectronNucleus(4,0) = 4.2500164887987364e-7;
    m_exactElectronNucleus(4,1) = 2.4927755071505254e-11;
    m_exactElectronNucleus(4,2) = 0.0049094042160749895;
    m_exactElectronNucleus(4,3) = -0.000502308838378577;
    m_exactElectronNucleus(4,4) = 0.0034477120078254845;
    m_exactElectronNucleus(4,5) = 1.5371142804860259e-24;
    m_exactElectronNucleus(4,6) = 4.924473721312616e-12;
    m_exactElectronNucleus(4,7) = 6.905163622343064e-20;
    m_exactElectronNucleus(4,8) = 8.30656406158157e-25;
    m_exactElectronNucleus(4,9) = 7.22099108679538e-14;
    m_exactElectronNucleus(5,0) = 0.0003964642966433137;
    m_exactElectronNucleus(5,1) = 0.003430101670771089;
    m_exactElectronNucleus(5,2) = 0.0008455744801229993;
    m_exactElectronNucleus(5,3) = 7.705198457410324e-6;
    m_exactElectronNucleus(5,4) = 2.4101085105907703e-24;
    m_exactElectronNucleus(5,5) = 0.00410220639105157;
    m_exactElectronNucleus(5,6) = 1.2630206441032246e-7;
    m_exactElectronNucleus(5,7) = 5.6014767862038486e-6;
    m_exactElectronNucleus(5,8) = 0.0008258291632664241;
    m_exactElectronNucleus(5,9) = -4.5500202831399375e-8;
    m_exactElectronNucleus(6,0) = 0.015936503246800314;
    m_exactElectronNucleus(6,1) = 0.0007350481705785905;
    m_exactElectronNucleus(6,2) = 0.009144848083587178;
    m_exactElectronNucleus(6,3) = 0.0328212556384292;
    m_exactElectronNucleus(6,4) = 1.3542259817811092e-12;
    m_exactElectronNucleus(6,5) = 7.016619072036424e-7;
    m_exactElectronNucleus(6,6) = 0.004711749527433483;
    m_exactElectronNucleus(6,7) = -4.7548817289418276e-6;
    m_exactElectronNucleus(6,8) = 3.1773271156945297e-9;
    m_exactElectronNucleus(6,9) = 4.6975115505053096e-5;
    m_exactElectronNucleus(7,0) = 0.000596978685186941;
    m_exactElectronNucleus(7,1) = -0.0012490246665673745;
    m_exactElectronNucleus(7,2) = 2.113669577326513e-5;
    m_exactElectronNucleus(7,3) = 7.57677376794402e-5;
    m_exactElectronNucleus(7,4) = 1.3633871200785802e-20;
    m_exactElectronNucleus(7,5) = -9.637059171293783e-6;
    m_exactElectronNucleus(7,6) = 8.967568581016109e-6;
    m_exactElectronNucleus(7,7) = 0.00020903847764887004;
    m_exactElectronNucleus(7,8) = -2.2328866337784158e-7;
    m_exactElectronNucleus(7,9) = 1.1471361708965405e-5;
    m_exactElectronNucleus(8,0) = 3.2917663138879515e-6;
    m_exactElectronNucleus(8,1) = 0.0003935800968818937;
    m_exactElectronNucleus(8,2) = -4.76367292089517e-5;
    m_exactElectronNucleus(8,3) = 1.3192670999483967e-6;
    m_exactElectronNucleus(8,4) = 2.0079402780434024e-24;
    m_exactElectronNucleus(8,5) = 0.00014134480259141638;
    m_exactElectronNucleus(8,6) = 7.0381882787075285e-9;
    m_exactElectronNucleus(8,7) = -9.807157201988673e-8;
    m_exactElectronNucleus(8,8) = 0.01717006262007232;
    m_exactElectronNucleus(8,9) = 1.3478816243453474e-10;
    m_exactElectronNucleus(9,0) = -3.910106601604258e-19;
    m_exactElectronNucleus(9,1) = -6.593276157701992e-5;
    m_exactElectronNucleus(9,2) = 0.000191339205906294;
    m_exactElectronNucleus(9,3) = 0.0005111756784408499;
    m_exactElectronNucleus(9,4) = 7.220991086795381e-14;
    m_exactElectronNucleus(9,5) = -8.051653540594755e-8;
    m_exactElectronNucleus(9,6) = 1.8301571728799462e-5;
    m_exactElectronNucleus(9,7) = 2.6909848608678254e-5;
    m_exactElectronNucleus(9,8) = 7.015232035588952e-12;
    m_exactElectronNucleus(9,9) = 0.00034822804501452684;
}

void IntegralTester::setupExactElectronElectronVector() {
    m_exactElectronElectron(0) = 0;
    m_exactElectronElectron(1) = 0.007166604041009602861;
    m_exactElectronElectron(2) = 0.022124581472837051566;
    m_exactElectronElectron(3) = 0.000138581030067768200;
    m_exactElectronElectron(4) = -6.8145328932903488e-08;
    m_exactElectronElectron(5) = 0.162484829269840000000;
    m_exactElectronElectron(6) = 0.266743478582807400000;
    m_exactElectronElectron(7) = 0.268120672073877200000;
    m_exactElectronElectron(8) = 0.526687299519774400000;
    m_exactElectronElectron(9) = -0.12730451834369380000;
}

IntegralTester::IntegralTester() {
    setupExactOverlapMatrix();
    setupExactKineticMatrix();
    setupExactElectronNucleusMatrix();
    setupExactElectronElectronVector();
    m_overlapIntegrator             = new OverlapIntegrator();
    m_kineticIntegrator             = new KineticIntegrator();
    m_electronNucleusIntegrator     = new ElectronNucleusIntegrator();
    m_electronElectronIntegrator    = new ElectronElectronIntegrator();

    m_primitives.reserve(10);
    m_primitives.push_back(GaussianPrimitive(0,0,0,     1.0,    vec{ 0.0, -1.5,  1.2}));
    m_primitives.push_back(GaussianPrimitive(1,0,0,     1.1,    vec{ 0.1, -1.2,  2.5}));
    m_primitives.push_back(GaussianPrimitive(0,1,0,     0.3,    vec{ 0.2, -1.0, -0.3}));
    m_primitives.push_back(GaussianPrimitive(0,0,1,     0.9,    vec{ 0.3, -0.8,  0.1}));
    m_primitives.push_back(GaussianPrimitive(0,0,2,     2.2,    vec{ 0.4, -0.6, -3.1}));
    m_primitives.push_back(GaussianPrimitive(0,2,0,     2.4,    vec{ 0.5, -0.4,  3.8}));
    m_primitives.push_back(GaussianPrimitive(2,0,0,     3.1,    vec{ 0.6, -0.2,  1.3}));
    m_primitives.push_back(GaussianPrimitive(1,1,0,     3.7,    vec{ 0.7,  0.0,  2.4}));
    m_primitives.push_back(GaussianPrimitive(0,1,1,     1.3,    vec{ 0.8,  0.2,  5.3}));
    m_primitives.push_back(GaussianPrimitive(1,0,1,     4.3,    vec{ 0.9,  0.4,  1.2}));
}

bool IntegralTester::runAllTests(bool silent) {
    bool passed = true;
    passed = runOverlapTests(silent)          && passed;
    passed = runKineticTests(silent)          && passed;
    passed = runElectronNucleusTests(silent)  && passed;
    passed = runElectronElectronTests(silent) && passed;

    if (passed) {
        if (!silent) cout << "   All integral tests PASSED." << endl;
    } else {
        if (!silent) cout << "   At least one integral test FAILED." << endl;
    }
    return passed;
}

bool IntegralTester::runOverlapTests(bool silent) {
    if (!silent) cout << "   Running overlap integral tests..." << endl;
    int failedTests = 0;
    int totalTests  = 100;

    for (int i=0; i<10; i++) {
        for (int j=0; j<10; j++) {
            double integral = m_overlapIntegrator->computeIntegral(&m_primitives[i], &m_primitives[j]);
            double difference = fabs(integral - m_exactOverlap(i,j));
            if (difference > m_tollerance) {
                failedTests++;
                if (!silent) cout << " Overlap integral test (" << i << "," << j << ") failed. Value: "
                                  << integral << ", exact value: " << m_exactOverlap(i,j) << ", abs. difference: "
                                  << difference << endl;
            }
        }
    }
    if (failedTests == 0) {
        if (!silent) cout << "      All 100 overlap integral tests PASSED." << endl;
    } else {
        if (!silent) cout << "      " << failedTests << " kinetic integral tests failed out of " << totalTests << "." << endl;
    }
    return (failedTests == 0 ? true : false);
}

bool IntegralTester::runKineticTests(bool silent) {
    if (!silent) cout << "   Running kinetic integral tests..." << endl;
    int failedTests = 0;
    int totalTests  = 100;

    for (int i=0; i<10; i++) {
        for (int j=0; j<10; j++) {
            double integral = m_kineticIntegrator->computeIntegral(&m_primitives[i], &m_primitives[j]);
            double difference = fabs(integral - m_exactKinetic(i,j));
            if (difference > m_tollerance) {
                failedTests++;
                if (!silent) cout << "      Kinetic integral test (" << i << "," << j << ") failed. Value: "
                                  << integral << ", exact value: " << m_exactKinetic(i,j) << ", abs. difference: "
                                  << difference << endl;
            }
        }
    }
    if (failedTests == 0) {
        if (!silent) cout << "      All 100 kinetic integral tests PASSED." << endl;
    } else {
        if (!silent) cout << "      " << failedTests << " kinetic integral tests failed out of " << totalTests << "." << endl;
    }
    return (failedTests == 0 ? true : false);
}

bool IntegralTester::runElectronNucleusTests(bool silent) {
    if (!silent) cout << "   Running electron-nucleus integral tests..." << endl;
    int failedTests = 0;
    int totalTests  = 100;

    for (int i=0; i<10; i++) {
        for (int j=0; j<10; j++) {
            int atom = (359*i+295*j-42*i*j+120*i*i-38*i*j*i) % 9;
            atom = (atom > 0 ? atom : -atom);
            vec nucleus = m_primitives[atom].nucleusPosition();
            double integral = m_electronNucleusIntegrator->computeIntegral(&m_primitives[i], &m_primitives[j], nucleus);
            double difference = fabs(integral - m_exactElectronNucleus(i,j));
            if (difference > m_tolleranceNumerical) {
                failedTests++;
                if (!silent) cout << "      Electron-nucleus integral test (" << i << "," << j << ") failed. Value: "
                                  << integral << ", exact value: " << m_exactElectronNucleus(i,j) << ", abs. difference: "
                                  << difference << endl;
            }
        }
    }
    if (failedTests == 0) {
        if (!silent) cout << "      All 100 electron-nucleus integral tests PASSED." << endl;
    } else {
        if (!silent) cout << "      " << failedTests << " electron-nucleus integral tests failed out of " << totalTests << "." << endl;
    }
    return (failedTests == 0 ? true : false);
}

bool IntegralTester::runElectronElectronTests(bool silent) {
    if (!silent) cout << "   Running electron-electron integral tests..." << endl;
    int failedTests = 0;
    int totalTests  = 9;

    for (int testIndex = 1; testIndex< 10; testIndex++) {
        double integral   = 0;
        double difference = 0;
        vec posA;
        vec posB;
        vec posC;
        vec posD;
        double a;
        double b;
        double c;
        double d;
        GaussianPrimitive primitiveA;
        GaussianPrimitive primitiveB;
        GaussianPrimitive primitiveC;
        GaussianPrimitive primitiveD;

        switch(testIndex) {
            case 1:
                posA = {-0.5, 0, 0};
                posB = {-0.5, 0, 0};
                posC = {-0.5, 0, 0};
                posD = {-0.5, 0, 0};
                a = 13.0077;
                b = 13.0077;
                c = 13.0077;
                d = 13.0077;
                primitiveA = GaussianPrimitive(0,0,0,a,posA);
                primitiveB = GaussianPrimitive(0,0,0,b,posB);
                primitiveC = GaussianPrimitive(0,0,0,c,posC);
                primitiveD = GaussianPrimitive(0,0,0,d,posD);
                break;

            case 2:
                posA = {0.5, 0, 0};
                posB = {-0.5, 0, 0};
                posC = {-0.5, 0, 0};
                posD = {0.5, 0, 0};
                a = 13.0077;
                b = 0.121949;
                c = 0.444529;
                d = 13.0077;
                primitiveA = GaussianPrimitive(0, 0, 0, a, posA);
                primitiveB = GaussianPrimitive(0, 0, 0, b, posB);
                primitiveC = GaussianPrimitive(0, 0, 0, c, posC);
                primitiveD = GaussianPrimitive(0, 0, 0, d, posD);
                break;

            case 3:
                posA = {0.5, 0, 0};
                posB = {-0.5, 0, 0};
                posC = {-0.5, 0, 0};
                posD = {0.5, 0, 0};
                a = 13.0077;
                b = 0.121949;
                c = 0.444529;
                d = 13.0077;
                primitiveA = GaussianPrimitive(0, 0 ,0, a, posA);
                primitiveB = GaussianPrimitive(0, 1 ,0, b, posB);
                primitiveC = GaussianPrimitive(0, 1 ,0, c, posC);
                primitiveD = GaussianPrimitive(0, 0 ,0, d, posD);
                break;

            case 4:
                posA = {0.55, 1, 3};
                posB = {-0.52, 5, 6};
                posC = {-0.53, 1, 2};
                posD = {0.45, 2, 4};
                a = 13.0077;
                b = 0.121949;
                c = 0.444529;
                d = 10.0077;
                primitiveA = GaussianPrimitive(1,0,0, a, posA);
                primitiveB = GaussianPrimitive(0,1,0, b, posB);
                primitiveC = GaussianPrimitive(0,1,0, c, posC);
                primitiveD = GaussianPrimitive(0,1,0, d, posD);
                break;

            case 5:
                posA = {1.2,2.3,3.4};
                posB = {-1.3,1.4,-2.4};
                posC = {2.3,0.9,3.2};
                posD = {5.0,1.9,1.2};
                a = 0.2;
                b = 0.3;
                c = 0.4;
                d = 0.1;
                primitiveA = GaussianPrimitive(0,0,0, a, posA);
                primitiveB = GaussianPrimitive(0,0,0, b, posB);
                primitiveC = GaussianPrimitive(0,0,0, c, posC);
                primitiveD = GaussianPrimitive(0,0,0, d, posD);
                break;

            case 6:
                posA = {1.2,2.3,3.4};
                posB = {-1.3,1.4,-2.4};
                posC = {2.3,0.9,3.2};
                posD = {5.0,1.9,1.2};
                a = 0.2;
                b = 0.3;
                c = 0.4;
                d = 0.1;
                primitiveA = GaussianPrimitive(0,0,0, a, posA);
                primitiveB = GaussianPrimitive(1,0,0, b, posB);
                primitiveC = GaussianPrimitive(0,0,0, c, posC);
                primitiveD = GaussianPrimitive(0,0,1, d, posD);
                break;

            case 7:
                posA = {1.2,2.3,3.4};
                posB = {-1.3,1.4,-2.4};
                posC = {2.3,0.9,3.2};
                posD = {5.0,1.9,1.2};
                a = 0.2;
                b = 0.3;
                c = 0.4;
                d = 0.1;
                primitiveA = GaussianPrimitive(0,0,0, a, posA);
                primitiveB = GaussianPrimitive(1,0,0, b, posB);
                primitiveC = GaussianPrimitive(0,2,0, c, posC);
                primitiveD = GaussianPrimitive(0,0,1, d, posD);
                break;

            case 8:
                posA = {1.2,2.3,3.4};
                posB = {-1.3,1.4,-2.4};
                posC = {2.3,0.9,3.2};
                posD = {5.0,1.9,1.2};
                a = 0.2;
                b = 0.3;
                c = 0.4;
                d = 0.1;
                primitiveA = GaussianPrimitive(1,1,0, a, posA);
                primitiveB = GaussianPrimitive(2,0,0, b, posB);
                primitiveC = GaussianPrimitive(2,0,0, c, posC);
                primitiveD = GaussianPrimitive(2,0,0, d, posD);
                break;

            case 9:
                posA = {1.2,2.3,3.4};
                posB = {-1.3,1.4,-2.4};
                posC = {2.3,0.9,3.2};
                posD = {5.0,1.9,1.2};
                a = 0.2;
                b = 0.3;
                c = 0.4;
                d = 0.1;
                primitiveA = GaussianPrimitive(1,1,0, a, posA);
                primitiveB = GaussianPrimitive(0,2,0, b, posB);
                primitiveC = GaussianPrimitive(0,2,0, c, posC);
                primitiveD = GaussianPrimitive(2,0,0, d, posD);
                break;

            default:
                break;
        }

        integral = m_electronElectronIntegrator->computeIntegral(&primitiveA,&primitiveB,&primitiveC,&primitiveD);
        difference = fabs(integral - m_exactElectronElectron(testIndex));
        if (difference > m_tollerance) {
            failedTests++;
            if (!silent) cout << "      Electron-nucleus integral test (" << testIndex << ") failed. Value: "
                              << integral << ", exact value: " << m_exactElectronElectron(testIndex) << ", abs. difference: "
                              << difference << endl;
        }

    }

    if (failedTests == 0) {
        if (!silent) cout << "      All 9 electron-electron integral tests PASSED." << endl;
    } else {
        if (!silent) cout << "      " << failedTests << " electron-electron integral tests failed out of " << totalTests << "." << endl;
    }
    return (failedTests == 0 ? true : false);
}
