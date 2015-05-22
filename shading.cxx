
struct LightingParameters
{
    LightingParameters(void)
    {
         lightDir[0] = +0.6;
         lightDir[1] = 0;
         lightDir[2] = +0.8;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 5.3;
         alpha = 7.5
    };
  

    doiuble lightDir[3]; // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};

LightingParameters lp;
