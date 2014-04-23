<?xml version="1.0" encoding="ISO-8859-1"?>

<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

    <xsl:output method="xml" indent="yes"/>

    <xsl:variable name="phaseOversampling">
        <xsl:choose>
            <xsl:when test="not(siemens/IRIS/DERIVED/phaseOversampling)">0</xsl:when>
            <xsl:otherwise>
                <xsl:value-of select="siemens/IRIS/DERIVED/phaseOversampling"/>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:variable>

    <xsl:variable name="sliceOversampling">
        <xsl:choose>
            <xsl:when test="not(siemens/MEAS/sKSpace/dSliceOversamplingForDialog)">0</xsl:when>
            <xsl:otherwise>
                <xsl:value-of select="siemens/MEAS/sKSpace/dSliceOversamplingForDialog"/>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:variable>

    <xsl:variable name="partialFourierPhase">
        <xsl:choose>
            <xsl:when test="siemens/MEAS/sKSpace/ucPhasePartialFourier = 1">0.5</xsl:when>
            <xsl:when test="siemens/MEAS/sKSpace/ucPhasePartialFourier = 2">0.75</xsl:when>
            <xsl:when test="siemens/MEAS/sKSpace/ucPhasePartialFourier = 4">0.875</xsl:when>
            <xsl:otherwise>1.0</xsl:otherwise>
        </xsl:choose>
    </xsl:variable>

    <xsl:variable name="numberOfContrasts">
        <xsl:value-of select="siemens/MEAS/lContrasts"/>
    </xsl:variable>


    <xsl:template match="/">
        <ismrmrdHeader xsi:schemaLocation="http://www.ismrm.org/ISMRMRD ismrmrd.xsd"
                xmlns="http://www.ismrm.org/ISMRMRD"
                xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
                xmlns:xs="http://www.w3.org/2001/XMLSchema">

            <!--
            <subjectInformation>
                <patientName>
                    <xsl:value-of select="siemens/HEADER/tPatientName"/>
                </patientName>
                <xsl:if test="siemens/YAPS/flUsedPatientWeight > 0">
                    <patientWeight_kg>
                        <xsl:value-of select="siemens/YAPS/flUsedPatientWeight"/>
                    </patientWeight_kg>
                </xsl:if>
            </subjectInformation>
            -->

            <acquisitionSystemInformation>
                <systemVendor>
                    <xsl:value-of select="siemens/DICOM/Manufacturer"/>
                </systemVendor>
                <systemModel>
                    <xsl:value-of select="siemens/DICOM/ManufacturersModelName"/>
                </systemModel>
                <systemFieldStrength_T>
                    <xsl:value-of select="siemens/YAPS/flMagneticFieldStrength"/>
                </systemFieldStrength_T>
                <receiverChannels>
                    <xsl:value-of select="siemens/YAPS/iMaxNoOfRxChannels" />
                </receiverChannels>
                <relativeReceiverNoiseBandwidth>0.79</relativeReceiverNoiseBandwidth>
            </acquisitionSystemInformation>

            <experimentalConditions>
                <H1resonanceFrequency_Hz>
                    <xsl:value-of select="siemens/DICOM/lFrequency"/>
                </H1resonanceFrequency_Hz>
            </experimentalConditions>

	    <!-- EPI Data Encoding Space 0 -->
            <encoding>
                <trajectory>epi</trajectory>
                <trajectoryDescription>
                  <identifier>ConventionalEPI</identifier>
                  <userParameterLong>
                    <name>etl</name>
                    <value>
                      <xsl:value-of select="siemens/MEAS/sFastImaging/lEPIFactor"/>
                    </value>
                  </userParameterLong>
                  <userParameterLong>
                    <name>numberOfNavigators</name>
                    <value>3</value>
                  </userParameterLong>
                  <userParameterLong>
                    <name>rampUpTime</name>
                    <value>
                      <xsl:value-of select="siemens/YAPS/alRegridRampupTime"/>
                    </value>
                  </userParameterLong>
                  <userParameterLong>
                    <name>rampDownTime</name>
                    <value>
                      <xsl:value-of select="siemens/YAPS/alRegridRampdownTime"/>
                    </value>
                  </userParameterLong>
                  <userParameterLong>
                    <name>flatTopTime</name>
                    <value>
                      <xsl:value-of select="siemens/YAPS/alRegridFlattopTime"/>
                    </value>
                  </userParameterLong>
                  <userParameterLong>
                    <name>echoSpacing</name>
                    <value>
                      <xsl:value-of select="siemens/YAPS/lEchoSpacing"/>
                    </value>
                  </userParameterLong>
                  <userParameterLong>
                    <name>acqDelayTime</name>
                    <value>
                      <xsl:value-of select="siemens/YAPS/alRegridDelaySamplesTime"/>
                    </value>
                  </userParameterLong>
                  <userParameterLong>
                    <name>numSamples</name>
                    <value>
                      <xsl:value-of select="siemens/YAPS/alRegridDestSamples"/>
                    </value>
                  </userParameterLong>
                  <userParameterDouble>
                    <name>dwellTime</name>
                    <value>
                      <xsl:value-of select="siemens/MEAS/sRXSPEC/alDwellTime div 1000.0"/>
                    </value>
                  </userParameterDouble>
                  <comment>Conventional 2D EPI sequence</comment>
                </trajectoryDescription>
                <encodedSpace>
                    <matrixSize>
		      <x>
			<xsl:value-of select="siemens/IRIS/DERIVED/imageColumns"/>
		      </x>
                      <y>
                        <xsl:value-of select="siemens/YAPS/iNoOfFourierLines"/>
                      </y>
                      <z>
                        <xsl:value-of select="siemens/YAPS/i3DFTLength"/>
                      </z>
                    </matrixSize>
                    <fieldOfView_mm>
		      <x>
			<xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dReadoutFOV"/>
                      </x>
                      <y>
                        <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dPhaseFOV * (1+$phaseOversampling)"/>
                      </y>
                      <z>
                        <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dThickness * (1+$sliceOversampling)"/>
                      </z>
                    </fieldOfView_mm>
                </encodedSpace>
                <reconSpace>
                    <matrixSize>
                        <x>
                            <xsl:value-of select="siemens/IRIS/DERIVED/imageColumns"/>
                        </x>
                        <y>
                            <xsl:value-of select="siemens/IRIS/DERIVED/imageLines"/>
                        </y>
                        <xsl:choose>
                            <xsl:when test="siemens/YAPS/i3DFTLength = 1">
                                <z>1</z>
                            </xsl:when>
                            <xsl:otherwise>
                                <z>
                                    <xsl:value-of select="siemens/MEAS/sKSpace/lImagesPerSlab"/>
                                </z>
                            </xsl:otherwise>
                        </xsl:choose>
                    </matrixSize>
                    <fieldOfView_mm>
                        <x>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dReadoutFOV"/>
                        </x>
                        <y>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dPhaseFOV"/>
                        </y>
                        <z>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dThickness"/>
                        </z>
                    </fieldOfView_mm>
                </reconSpace>
                <encodingLimits>
                    <kspace_encoding_step_1>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:value-of select="siemens/YAPS/iNoOfFourierLines - 1"/>
                        </maximum>
                        <xsl:choose>
                            <xsl:when test="siemens/MEAS/sKSpace/ucTrajectory = 1">
                                <center>
                                    <xsl:value-of select="floor((siemens/MEAS/sKSpace/lPhaseEncodingLines div 2)-(siemens/MEAS/sKSpace/lPhaseEncodingLines - siemens/YAPS/iNoOfFourierLines))"/>
                                </center>
                            </xsl:when>
                            <xsl:otherwise>
                                <center>0</center>
                            </xsl:otherwise>
                        </xsl:choose>
                    </kspace_encoding_step_1>
                    <kspace_encoding_step_2>
                        <minimum>0</minimum>
                        <xsl:choose>
                            <xsl:when test="not(siemens/YAPS/iNoOfFourierPartitions) or (siemens/YAPS/i3DFTLength = 1)">
                                <maximum>0</maximum>
                                <center>0</center>
                            </xsl:when>
                            <xsl:otherwise>
                                <maximum>
                                    <xsl:value-of select="siemens/YAPS/iNoOfFourierPartitions - 1"/>
                                </maximum>
                                <center>
                                    <xsl:value-of select="floor((siemens/MEAS/sKSpace/lPartitions div 2)-(siemens/MEAS/sKSpace/lPartitions - siemens/YAPS/iNoOfFourierPartitions))"/>
                                </center>
                            </xsl:otherwise>
                        </xsl:choose>
                    </kspace_encoding_step_2>
                    <slice>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/lSize - 1"/>
                        </maximum>
                        <center>0</center>
                    </slice>
                    <repetition>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:choose>
                                <xsl:when test="siemens/MEAS/lRepetitions">
                                    <xsl:value-of select="siemens/MEAS/lRepetitions - 1"/>
                                </xsl:when>
                                <xsl:otherwise>0</xsl:otherwise>
                            </xsl:choose>
                        </maximum>
                        <center>0</center>
                    </repetition>
                    <contrast>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:choose>
                                <xsl:when test="siemens/MEAS/lContrasts">
                                    <xsl:value-of select="siemens/MEAS/lContrasts - 1"/>
                                </xsl:when>
                                <xsl:otherwise>0</xsl:otherwise>
                            </xsl:choose>
                        </maximum>
                        <center>0</center>
                    </contrast>
                </encodingLimits>
            </encoding>

	    <!-- FLASH PAT Reference Scan Encoding Space 1-->
            <encoding>
                <trajectory>cartesian</trajectory>
                <encodedSpace>
                    <matrixSize>
		      <x><xsl:value-of select="siemens/IRIS/DERIVED/imageColumns"/></x>
                      <y><xsl:value-of select="siemens/YAPS/iNoOfFourierLines"/></y>
                      <z><xsl:value-of select="siemens/YAPS/i3DFTLength"/></z>
                    </matrixSize>
                    <fieldOfView_mm>
		      <x><xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dReadoutFOV"/></x>
                      <y><xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dPhaseFOV * (1+$phaseOversampling)"/></y>
                      <z><xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dThickness * (1+$sliceOversampling)"/></z>
                    </fieldOfView_mm>
                </encodedSpace>
                <reconSpace>
                    <matrixSize>
                        <x><xsl:value-of select="siemens/IRIS/DERIVED/imageColumns"/></x>
                        <y><xsl:value-of select="siemens/IRIS/DERIVED/imageLines"/></y>
                        <xsl:choose>
                            <xsl:when test="siemens/YAPS/i3DFTLength = 1">
                                <z>1</z>
                            </xsl:when>
                            <xsl:otherwise>
                                <z><xsl:value-of select="siemens/MEAS/sKSpace/lImagesPerSlab"/></z>
                            </xsl:otherwise>
                        </xsl:choose>
                    </matrixSize>
                    <fieldOfView_mm>
                        <x><xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dReadoutFOV"/></x>
                        <y><xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dPhaseFOV"/></y>
                        <z><xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dThickness"/></z>
                    </fieldOfView_mm>
                </reconSpace>
                <encodingLimits>
                    <kspace_encoding_step_1>
                        <minimum>0</minimum>
                        <maximum><xsl:value-of select="siemens/YAPS/iNoOfFourierLines - 1"/></maximum>
                        <center>
			  <xsl:value-of select="floor((siemens/MEAS/sKSpace/lPhaseEncodingLines div 2)-(siemens/MEAS/sKSpace/lPhaseEncodingLines - siemens/YAPS/iNoOfFourierLines))"/>
                        </center>
                    </kspace_encoding_step_1>
                    <kspace_encoding_step_2>
                        <minimum>0</minimum>
                        <xsl:choose>
                            <xsl:when test="not(siemens/YAPS/iNoOfFourierPartitions) or (siemens/YAPS/i3DFTLength = 1)">
                                <maximum>0</maximum>
                                <center>0</center>
                            </xsl:when>
                            <xsl:otherwise>
                                <maximum>
                                    <xsl:value-of select="siemens/YAPS/iNoOfFourierPartitions - 1"/>
                                </maximum>
                                <center>
                                    <xsl:value-of select="floor((siemens/MEAS/sKSpace/lPartitions div 2)-(siemens/MEAS/sKSpace/lPartitions - siemens/YAPS/iNoOfFourierPartitions))"/>
                                </center>
                            </xsl:otherwise>
                        </xsl:choose>
                    </kspace_encoding_step_2>
                    <slice>
                        <minimum>0</minimum>
                        <maximum><xsl:value-of select="siemens/MEAS/sSliceArray/lSize - 1"/></maximum>
                        <center>0</center>
                    </slice>
                </encodingLimits>
            </encoding>

            <parallelImaging>
                <accelerationFactor>
                    <kspace_encoding_step_1>
                        <xsl:choose>
                            <xsl:when test="not(siemens/MEAS/sPat/lAccelFactPE)">1</xsl:when>
                            <xsl:otherwise>
                                <xsl:value-of select="(siemens/MEAS/sPat/lAccelFactPE)"/>
                            </xsl:otherwise>
                        </xsl:choose>
                    </kspace_encoding_step_1>
                    <kspace_encoding_step_2>
                        <xsl:choose>
                            <xsl:when test="not(siemens/MEAS/sPat/lAccelFact3D)">1</xsl:when>
                            <xsl:otherwise>
                                <xsl:value-of select="(siemens/MEAS/sPat/lAccelFact3D)"/>
                            </xsl:otherwise>
                        </xsl:choose>
                    </kspace_encoding_step_2>
                </accelerationFactor>
                <calibrationMode>
                    <xsl:choose>
                        <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 1">other</xsl:when>
                        <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 2">embedded</xsl:when>
                        <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 4">separate</xsl:when>
                        <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 8">separate</xsl:when>
                        <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 16">interleaved</xsl:when>
                        <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 32">interleaved</xsl:when>
                        <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 64">interleaved</xsl:when>
                        <xsl:otherwise>other</xsl:otherwise>
                    </xsl:choose>
                </calibrationMode>
                <xsl:if test="(siemens/MEAS/sPat/ucRefScanMode = 1) or (siemens/MEAS/sPat/ucRefScanMode = 16) or (siemens/MEAS/sPat/ucRefScanMode = 32) or (siemens/MEAS/sPat/ucRefScanMode = 64)">
                    <interleavingDimension>
                        <xsl:choose>
                            <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 16">average</xsl:when>
                            <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 32">repetition</xsl:when>
                            <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 64">phase</xsl:when>
                            <xsl:otherwise>other</xsl:otherwise>
                        </xsl:choose>
                    </interleavingDimension>
                </xsl:if>
            </parallelImaging>

            <sequenceParameters>
                <xsl:for-each select="siemens/MEAS/alTR">
                    <xsl:if test=". &gt; 0">
                        <TR>
                            <xsl:value-of select=". div 1000.0" />
                        </TR>
                    </xsl:if>
                </xsl:for-each>
                <xsl:for-each select="siemens/MEAS/alTE">
                    <xsl:if test=". &gt; 0">
                        <xsl:if test="position() &lt; ($numberOfContrasts + 1)">
                            <TE>
                                <xsl:value-of select=". div 1000.0" />
                            </TE>
                        </xsl:if>
                    </xsl:if>
                </xsl:for-each>
                <xsl:for-each select="siemens/MEAS/alTI">
                    <xsl:if test=". &gt; 0">
                        <TI>
                            <xsl:value-of select=". div 1000.0" />
                        </TI>
                    </xsl:if>
                </xsl:for-each>
            </sequenceParameters>

        </ismrmrdHeader>
    </xsl:template>

</xsl:stylesheet>
