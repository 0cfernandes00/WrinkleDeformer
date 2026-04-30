__kernel void wrinkleDeform(
    __global float* finalPositions,
    __global const float* initialPositions,
    __global const float* vertexNormals,  // needs to be added
    __global const float* wrinklePhase,   // uploaded from CPU BFS result
    __global const float* strainMask,
    __global const float* vertexAmps,
    const float envelope,
    const float compressionThreshold,
    const uint numElements)
{
    uint id = get_global_id(0);
    if (id >= numElements) return;
    uint p = id * 3;

    float strain = strainMask[id];
    if (strain > 0.0001f) {
        float phase = wrinklePhase[id];
        float amp = vertexAmps[id];
        float wave = pow(fabs(sin(phase)), 0.5f) * (sin(phase) > 0.0f ? 1.0f : -1.0f);
        float disp = wave * amp * envelope;

        finalPositions[p] = initialPositions[p] + vertexNormals[p] * disp;
        finalPositions[p + 1] = initialPositions[p + 1] + vertexNormals[p + 1] * disp;
        finalPositions[p + 2] = initialPositions[p + 2] + vertexNormals[p + 2] * disp;
    }
    else {
        finalPositions[p] = initialPositions[p];
        finalPositions[p + 1] = initialPositions[p + 1];
        finalPositions[p + 2] = initialPositions[p + 2];
    }
}