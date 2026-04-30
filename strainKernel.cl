__kernel void wrinkleStrain(
	__global const float* triQInv,           // [numTris * 4] — qi00,qi01,qi10,qi11
	__global const float* triNormals,        // [numTris * 6] — n1.xyz, n2.xyz
	__global const int* triVertIdx,          // [numTris * 3]
	__global const float* triRestCrossSqLen, // [numTris]
	
	__global float* vertStrainAccum,         // [numVerts * 3] — x=strainMag, y=physAmp, z=count
	__global float* vertDirAccum,            // [numVerts * 3] — wrinkle dir xyz
	__global float* vertAreaDelta           // [numVerts * 3]

) {

}
