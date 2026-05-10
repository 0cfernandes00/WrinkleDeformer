void atomic_add_float(__global volatile float* addr, float val) {
    union { float f; uint i; } old, newval;
    do {
        old.f = *addr;
        newval.f = old.f + val;
    } while (atom_cmpxchg((__global volatile uint*)addr, old.i, newval.i) != old.i);
}

__kernel void wrinkleStrain(
    __global const float* initialPositions,  // [numVerts * 3]
    __global const float* triQInv,            // [numTris * 4]  qi00,qi01,qi10,qi11
    __global const float* triNormals,         // [numTris * 6]  n1.xyz, n2.xyz
    __global const int* triVertIdx,         // [numTris * 3]
    __global const float* triRestCrossSqLen,  // [numTris]
    __global volatile float* vertStrainAccum,    // [numVerts * 3]  strainMag, physAmp, count
    __global volatile float* vertDirAccum,       // [numVerts * 3]  wrinkle dir xyz
    __global volatile float* vertAreaDelta,      // [numVerts * 3]
    const float warpStiffness,
    const float weftStiffness,
    const float areaStiffness,
    const float frequency,
    const float amplitude,
    const uint  numTris
)
{
    uint id = get_global_id(0);
    if (id >= numTris) return;

    // -------------------------------------------------------------------------
    // Unpack triangle topology
    // -------------------------------------------------------------------------
    float qi00 = triQInv[id * 4 + 0], qi01 = triQInv[id * 4 + 1];
    float qi10 = triQInv[id * 4 + 2], qi11 = triQInv[id * 4 + 3];

    int v0 = triVertIdx[id * 3 + 0];
    int v1 = triVertIdx[id * 3 + 1];
    int v2 = triVertIdx[id * 3 + 2];

    float n1x = triNormals[id * 6 + 0], n1y = triNormals[id * 6 + 1], n1z = triNormals[id * 6 + 2];
    float n2x = triNormals[id * 6 + 3], n2y = triNormals[id * 6 + 4], n2z = triNormals[id * 6 + 5];

    float restCrossSqLen = triRestCrossSqLen[id];

    // -------------------------------------------------------------------------
    // Read current vertex positions
    // -------------------------------------------------------------------------
    float ax = initialPositions[v0 * 3], ay = initialPositions[v0 * 3 + 1], az = initialPositions[v0 * 3 + 2];
    float bx = initialPositions[v1 * 3], by = initialPositions[v1 * 3 + 1], bz = initialPositions[v1 * 3 + 2];
    float cx = initialPositions[v2 * 3], cy = initialPositions[v2 * 3 + 1], cz = initialPositions[v2 * 3 + 2];

    // -------------------------------------------------------------------------
    // P matrix — columns are (p1-p0) and (p2-p0)
    // -------------------------------------------------------------------------
    float p00 = bx - ax, p01 = cx - ax;
    float p10 = by - ay, p11 = cy - ay;
    float p20 = bz - az, p21 = cz - az;

    // -------------------------------------------------------------------------
    // F = P * Q^-1  (deformation gradient)
    // -------------------------------------------------------------------------
    float f00 = p00 * qi00 + p01 * qi10;
    float f01 = p00 * qi01 + p01 * qi11;
    float f10 = p10 * qi00 + p11 * qi10;
    float f11 = p10 * qi01 + p11 * qi11;
    float f20 = p20 * qi00 + p21 * qi10;
    float f21 = p20 * qi01 + p21 * qi11;

    float f1Len = sqrt(f00 * f00 + f10 * f10 + f20 * f20);
    float f2Len = sqrt(f01 * f01 + f11 * f11 + f21 * f21);

    // -------------------------------------------------------------------------
    // S = F^T * F  (right Cauchy-Green / strain tensor, symmetric 2x2)
    // -------------------------------------------------------------------------
    float s00 = f00 * f00 + f10 * f10 + f20 * f20;
    float s01 = f00 * f01 + f10 * f11 + f20 * f21;
    float s11 = f01 * f01 + f11 * f11 + f21 * f21;
    // s10 == s01 by symmetry

    float s01_decoupled = (f1Len * f2Len > 1e-6f) ? s01 / (f1Len * f2Len) : 0.f;

    float strainMagnitude = (s00 - 1.0f) * warpStiffness
        + (s11 - 1.0f) * weftStiffness;

    // -------------------------------------------------------------------------
    // Area constraint  C_area = |e1 x e2|^2 - restCrossSqLen
    // e1 = col0 of P = (p00,p10,p20),  e2 = col1 = (p01,p11,p21)
    // -------------------------------------------------------------------------
    float cross_x = p10 * p21 - p20 * p11;
    float cross_y = p20 * p01 - p00 * p21;
    float cross_z = p00 * p11 - p10 * p01;
    float crossSqLen = cross_x * cross_x + cross_y * cross_y + cross_z * cross_z;
    float C_area = crossSqLen - restCrossSqLen;

    if (C_area < 0.0f) {
        // grad1 = 2 * e2 x cross_p
        float g1x = 2.0f * (p11 * cross_z - p21 * cross_y);
        float g1y = 2.0f * (p21 * cross_x - p01 * cross_z);
        float g1z = 2.0f * (p01 * cross_y - p11 * cross_x);

        // grad2 = 2 * e1 x (e2 x e1) = 2 * e1 x (-cross_p)
        float g2x = 2.0f * (-p10 * cross_z + p20 * cross_y);
        float g2y = 2.0f * (-p20 * cross_x + p00 * cross_z);
        float g2z = 2.0f * (-p00 * cross_y + p10 * cross_x);

        float g0x = -(g1x + g2x);
        float g0y = -(g1y + g2y);
        float g0z = -(g1z + g2z);

        float denom = g0x * g0x + g0y * g0y + g0z * g0z
            + g1x * g1x + g1y * g1y + g1z * g1z
            + g2x * g2x + g2y * g2y + g2z * g2z;

        if (denom >= 1e-8f) {
            float s = -(C_area / denom) * areaStiffness;

            atomic_add_float(&vertAreaDelta[v0 * 3 + 0], g0x * s);
            atomic_add_float(&vertAreaDelta[v0 * 3 + 1], g0y * s);
            atomic_add_float(&vertAreaDelta[v0 * 3 + 2], g0z * s);
            atomic_add_float(&vertAreaDelta[v1 * 3 + 0], g1x * s);
            atomic_add_float(&vertAreaDelta[v1 * 3 + 1], g1y * s);
            atomic_add_float(&vertAreaDelta[v1 * 3 + 2], g1z * s);
            atomic_add_float(&vertAreaDelta[v2 * 3 + 0], g2x * s);
            atomic_add_float(&vertAreaDelta[v2 * 3 + 1], g2y * s);
            atomic_add_float(&vertAreaDelta[v2 * 3 + 2], g2z * s);
        }
    }

    // -------------------------------------------------------------------------
    // 2x2 eigenvalue decomposition of S for compression factor
    // -------------------------------------------------------------------------
    float trace = s00 + s11;
    float det = s00 * s11 - s01_decoupled * s01_decoupled;
    if (fabs(det) < 1e-6f) return;

    float disc = sqrt(fmax(0.0f, (trace * trace * 0.25f) - det));
    float eigenMin = (trace * 0.5f) - disc;
    float compFactor = sqrt(fmax(0.0f, eigenMin));

    // -------------------------------------------------------------------------
    // Physical wrinkle amplitude
    // -------------------------------------------------------------------------
    float physAmp = 0.0f;
    if (compFactor >= 0.0001f && compFactor < 1.0f) {
        float waveLen = 1.0f / frequency;
        float invC = 1.0f / compFactor;
        physAmp = (waveLen / (2.0f * 3.14159265358979f))
            * sqrt(fmax(0.0f, invC * invC - 1.0f));
    }
    physAmp *= amplitude;

    // -------------------------------------------------------------------------
    // Wrinkle direction (minimum principal strain axis → perpendicular in UV)
    // -------------------------------------------------------------------------
    float angle = 0.5f * atan2(2.0f * s01_decoupled, s00 - s11);
    float compX = cos(angle);
    float compY = sin(angle);
    float wrinkleX = -compY;
    float wrinkleY = compX;

    // Transform UV-space direction to world space via stored tangent frame
    float wdx = n1x * wrinkleX + n2x * wrinkleY;
    float wdy = n1y * wrinkleX + n2y * wrinkleY;
    float wdz = n1z * wrinkleX + n2z * wrinkleY;
    float dirLen = sqrt(wdx * wdx + wdy * wdy + wdz * wdz);
    if (dirLen < 1e-6f) return;
    wdx /= dirLen;
    wdy /= dirLen;
    wdz /= dirLen;

    // -------------------------------------------------------------------------
    // Scatter to all 3 vertices via atomics
    // vertStrainAccum layout: [strainMag, physAmp, count] per vertex
    // Direction flip (coherence check) is omitted — the atomic scatter pattern
    // makes read-modify-write coherence unsafe; averaging converges well enough.
    // -------------------------------------------------------------------------
    int verts[3];
    verts[0] = v0; verts[1] = v1; verts[2] = v2;

    for (int j = 0; j < 3; j++) {
        int v = verts[j];
        atomic_add_float(&vertStrainAccum[v * 3 + 0], strainMagnitude);
        atomic_add_float(&vertStrainAccum[v * 3 + 1], physAmp);
        atomic_add_float(&vertStrainAccum[v * 3 + 2], 1.0f);  // count for averaging
        atomic_add_float(&vertDirAccum[v * 3 + 0], wdx);
        atomic_add_float(&vertDirAccum[v * 3 + 1], wdy);
        atomic_add_float(&vertDirAccum[v * 3 + 2], wdz);
    }
}
