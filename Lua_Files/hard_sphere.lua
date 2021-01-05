simion.workbench_program()

-- The code implements a hard-sphere collision model.

-- Features and assumptions of the model:
-- - Ion collisions follow the hard-sphere collision model.
--     Energy transfers occur solely via these collisions.
-- - Ion collisions are elastic.
-- - Background gas is assumed neutral in charge.
-- - Background gas velocity follows the Maxwell-Boltzmann distribution.
-- - Background gas mean velocity may be non-zero.
-- - Kinetic cooling and heating of ions due to collisions are simulated.
-- - Kinetic cooling and heating of background gas is assumed negligible over many collisions.
--
-- Note on time-steps: each individual ion-gas collision is modeled, which requires the time-step to be some fraction of mean-free-path.
-- Therefore, simulations with frequent collisions (i.e. higher pressure) can be computationally intensive.

io.input("C:\\Users\\aydan\\SIMION\\Garcia_Lab\\PAs\\current_params.txt")

adjustable ac_voltage = tonumber(io.read("*line"))   -- AC voltage
adjustable dc_voltage = tonumber(io.read("*line"))   -- DC voltage
adjustable omega = 1000 * 3.14      -- angular frequency (radians per microsecond)

-- Mean free path (MFP) (mm) between collisions.
-- Set to -1 (recommended) to calculate this automatically from pressure and temperature.
adjustable _mean_free_path_mm = 10

-- Mass of background gas particle (amu)
adjustable _gas_mass_amu = 2
adjustable _temperature_k = 300
adjustable _pressure_pa = 0.53

-- Collision-cross section (m^2)
-- (The diameter of the cross-sectional area is roughly
--  the sum of the diameters of the colliding ion and buffer gas particles.)
-- (2.1E-19 is roughly for two Helium atoms--Atkins1998-Table 1.3)
-- (Note: the Van der Waals radius for He is 140 pm = 1.40 angstrom.
--   -- http://www.webelements.com and http://en.wikipedia.org/wiki/Helium --
--   i.e. 2.46e-19 collision cross-section)
-- (2.27E-18 is for collision between He and some 200 amu ion with combined
--  collision diameter of 2 + 15 angstroms.  It is used in some benchmarks.)
adjustable _sigma_m2 = 2.27E-18

-- Mean background gas velocity (mm/usec) in x,y,z directions.
-- Normally, these are zero.
adjustable _vx_bar_gas_mmusec = 0
adjustable _vy_bar_gas_mmusec = 0
adjustable _vz_bar_gas_mmusec = 0

-- Mean number of time steps per MFP.
-- Typically this default is ok.  We want sufficient number of
-- time-steps per mean-free path for this code to be reliable.
adjustable _steps_per_MFP = 20.0

-- Collision marker flag.
-- If non-zero, markers will be placed at the collisions.
adjustable _mark_collisions = 1

-- How much trace data (average KE) to output.
-- (0=none, 1=at each splat, 2=at each collision)
adjustable _trace_level = 0
-- If _trace_level is 2, this is the number of collisions before each trace.
-- This reduces the verbosity of the trace.
adjustable _trace_skip = 100

---- Internal variables

-- Statistics

---- current running average of KE for each particle.  maps ion_number --> KE.
local ke_averages = {}

---- last collision time for each particle.  maps ion_number --> time.
local last_collision_times = {}

-- Last known ion number (-1 = undefined).
local last_ion_number = -1

-- Last known ion speed (-1 = undefined).
local last_speed_ion = -1

-- Currently used mean-free path (-1 = undefined).
local effective_mean_free_path_mm = -1

-- Count relative to _trace_skip
local trace_count = 0

-- Maximum time step (usec) that fast_adjust should permit.
-- This is continually updated so that the _steps_per_MFP setting
-- remains meaningful.
local max_timestep

-- Define constants
local k = 1.3806505e-23       -- Boltzmann constant (J/K)
-- local R = 8.3145           -- Ideal gas constant (J/(mol*K))
local kg_amu = 1.6605402e-27  -- (kg/amu) conversion factor
local pi = 3.1415926535       -- PI constant
local eV_J = 6.2415095e+18    -- (eV/J) conversion factor


-- Error function (erf).
--   erf(z) = (2/sqrt(pi)) * integral[0..z] exp(-t^2) dt
-- This algorithm is quite accurate.  It is based on
-- "Q-Function Handout" by Keith Chugg:
--   http://tesla.csl.uiuc.edu/~koetter/ece361/Q-function.pdf
-- See also http://www.theorie.physik.uni-muenchen.de/~serge/erf-approx.pdf
-- I also find that the following makes a reasonable approximation:
--   1 - exp(-(2/sqrt(pi))x - (2/pi)x^2)
function erf(z)
    local z2 = abs(z)
    local t = 1 / (1 + 0.32759109962 * z2)
    local res = (    - 1.061405429 ) * t
    res = (res + 1.453152027 ) * t
    res = (res - 1.421413741 ) * t
    res = (res + 0.2844966736) * t
    res =((res - 0.254829592 ) * t) * exp(-z2*z2)
    res = res + 1
    if z < 0 then res = -res end
    return res
end

-- Return a normalized Gaussian random variable (-inf, +inf).
-- [ http://en.wikipedia.org/wiki/Normal_distribution ]
function gaussian_random()
    -- Using the Box-Muller algorithm.
    local s = 1
    local v1, v2
    while s >= 1 do
        v1 = 2*rand() - 1
        v2 = 2*rand() - 1
        s = v1*v1 + v2*v2
    end
    local rand1 = v1*sqrt(-2*ln(s) / s)  -- (assume divide by zero improbable?)
    return rand1
end

-- SIMION initialize segment. Called on particle creation.
function segment.initialize()
end

-- SIMION time step adjust segment. Adjusts time step sizes.
function segment.tstep_adjust()
    -- Ensure time steps are sufficiently small.  They should be some
    -- fraction of mean-free-path so that collisions are not missed.
    if max_timestep and ion_time_step > max_timestep then
        ion_time_step = max_timestep
    end
end

-- SIMION other actions segment. Called on every time step.
function segment.other_actions()
    if _pressure_pa == 0 then  -- collisions disabled
        return
    end

    -- Temporarily translate ion velocity (mm/us) frame of
    -- reference such that mean background gas velocity is zero.
    -- This simplifies the subsequent analysis.
    local vx = ion_vx_mm - _vx_bar_gas_mmusec
    local vy = ion_vy_mm - _vy_bar_gas_mmusec
    local vz = ion_vz_mm - _vz_bar_gas_mmusec

    -- Obtain ion speed (relative to mean background gas velocity).
    local speed_ion = sqrt(vx^2 + vy^2 + vz^2)
    if speed_ion < 1E-7 then
         speed_ion = 1E-7  -- prevent divide by zero and such effects later on
    end

    -- Compute mean-free-path.
    -- > See notes.pdf for discussion on the math.
    if _mean_free_path_mm > 0 then -- explicitly specified
        effective_mean_free_path_mm = _mean_free_path_mm
    else  -- calculate from current ion velocity
        -- Only recompute mean-free-path if speed_ion has
        -- changed significantly. This is intended to speed up the
        -- calculation a bit.  This code will handle flying ions by groups.
        if last_ion_number ~= ion_number or
                abs(speed_ion / last_speed_ion - 1) > 0.05  -- changed
        then
            -- Compute mean gas speed (mm/us)
            local c_bar_gas = sqrt(8*k*_temperature_k/pi/(_gas_mass_amu * kg_amu)) / 1000

            -- Compute median gas speed (mm/us)
            local c_star_gas = sqrt(2*k*_temperature_k/(_gas_mass_amu * kg_amu)) / 1000

            -- Compute mean relative speed (mm/us) between gas and ion.
            local s = speed_ion / c_star_gas
            local c_bar_rel = c_bar_gas * (
                (s + 1/(2*s)) * 0.5 * sqrt(pi) * erf(s) + 0.5 * exp(-s*s))

            -- Compute mean-free-path (mm)
            effective_mean_free_path_mm = 1000 * k * _temperature_k *
                (speed_ion / c_bar_rel) / (_pressure_pa * _sigma_m2)

            -- Store data about this calculation.
            last_speed_ion = speed_ion
            last_ion_number = ion_number

            --print("DEBUG:ion[c],gas[c_bar],c_bar_rel,MFP=",
            --      speed_ion, c_bar_gas, c_bar_rel, effective_mean_free_path_mm)

            -- Note: The following is a simpler and almost as suitable
            -- approximation for c_bar_rel, which you may used instead:
            -- c_bar_rel = sqrt(speed_ion^2 + c_bar_gas^2)
        end
    end

    -- Limit time-step size to a fraction of the MFP.
    max_timestep = effective_mean_free_path_mm / speed_ion / _steps_per_MFP

    -- Compute probability of collision in current time-step.
    -- > For an infinitesimal distance (dx) traveled, the increase in the
    --   fraction (f) of collided particles relative to the number
    --   of yet uncollided particles (1-f) is equal to the distance
    --   traveled (dx) over the mean-free-path (lambda):
    --     df/(1-f) = dx / lambda
    --   Solving this differential equation gives
    --     f = 1 - exp(- dx / lambda) = 1 - exp(- v dt / lambda)
    --   This f can be interpreted as the probability that a single
    --   particle collides in the distance traveled.
    local collision_prob = 1 -
        exp(- speed_ion * ion_time_step / effective_mean_free_path_mm)

    -- Test for collision.
    if rand() > collision_prob then
        return -- no collision
    end

    ----- Handle collision.

    -- Compute standard deviation of background gas velocity in
    -- one dimension (mm/us).
    -- > From kinetic gas theory (Maxwell-Boltzmann), velocity in
    --   one dimension is normally distributed with standard
    --   deviation sqrt(kT/m).
    local vr_stdev_gas =
        sqrt(k * _temperature_k / (_gas_mass_amu * kg_amu)) / 1000

    -- Compute velocity of colliding background gas particle.
    -- > For the population of background gas particles that collide with the
    --   ion, their velocities are not entirely Maxwell (Gaussian) but
    --   are also proportional to the relative velocities the ion and
    --   background gas particles:
    --     p(v_gas) = |v_gas - v_ion| f(v_gas)
    --   See notes.pdf for discussion.
    -- > To generate random velocities in this distribution, we may
    --   use a rejection method (http://en.wikipedia.org/wiki/Rejection_sampling)
    --   approach:
    --   > Pick a gas velocity from the Maxwell distribution.
    --   > Accept with probability proportional to its
    --     speed relative to the ion.
    local vx_gas, vy_gas, vz_gas -- computed velocities
    -- > scale is an approximate upper-bound for "len" calculated below.
    --   We'll use three standard deviations of the three dimensional gas velocity.
    local scale = speed_ion + vr_stdev_gas * 1.732 * 3  --sqrt(3)=~1.732
    repeat
        vx_gas = gaussian_random() * vr_stdev_gas
        vy_gas = gaussian_random() * vr_stdev_gas
        vz_gas = gaussian_random() * vr_stdev_gas
        local len = sqrt((vx_gas - vx)^2 + (vy_gas - vy)^2 + (vz_gas - vz)^2)
        --assert(len <= scale) -- true at least ~99% of the time.
    until rand() < len / scale

    -- Alernately, for greater performance and as an approximation, you might
    -- replace the above with a simple Maxwell distribution:
    --  vx_gas = gaussian_random() * vr_stdev_gas
    --  vy_gas = gaussian_random() * vr_stdev_gas
    --  vz_gas = gaussian_random() * vr_stdev_gas

    -- Translate velocity reference frame so that colliding
    -- background gas particle is stationary.
    -- > This simplifies the subsequent analysis.
    vx = vx - vx_gas
    vy = vy - vy_gas
    vz = vz - vz_gas

    -- > Notes on collision orientation
    --   A collision of the ion in 3D can now be reasoned in 2D since
    --   the ion remains in some 2D plane before and after collision.
    --   The ion collides with an gas particle initially at rest (in the
    --   current velocity reference frame).
    --   For convenience, we define a coordinate system (r, t) on the
    --   collision plane.  r is the radial axis through the centers of
    --   the colliding particles, with the positive direction indicating
    --   approaching particles.  t is the tangential axis perpendicular to r.
    --   An additional coordinate theta defines the the rotation of the
    --   collision plane around the ion velocity axis.

    -- Compute randomized impact offset [0, 1) as a fraction
    -- of collisional cross-section diameter.
    -- 0 is a head-on collision; 1 would be a near miss.
    -- > You can imaging this as the gas particle being a stationary
    --   dart board of radius 1 unit (representing twice the actual radius
    --   of the gas particle) and the ion center is a dart
    --   with velocity perpendicular to the dart board.
    --   The dart has equal probability of hitting anywhere on the
    --   dart board.  Since a radius "d" from the center represents
    --   a ring with circumference proportional to "d", this implies
    --   that the probability of hitting at a distance "d" from the
    --   center is proportional to "d".
    -- > Formally, the normalized probability density function is
    --   f(d) = 2*d for d in [0,1].  From the fundamental transformation
    --   law of probabilities, we have
    --   integral[0..impact_offset] f(d) dd = impact_offset^2 = U,
    --   where U is a uniform random variable.  That is,
    --   impact_offset = sqrt(U).  Decrease it it slightly
    --   to prevent overflow in asin later.
    local impact_offset = sqrt(0.999999999 * rand())

    -- Convert impact offset to impact angle [0, +pi/2) (radians).
    -- Do this since the target is really a sphere (not flat dartboard).
    -- This is the angle between the relative velocity
    -- between the two colliding particles (i.e. the velocity of the dart
    -- imagined perpendicular to the dart board) and the r axis
    -- (i.e. a vector from the center of the gas particle to the location
    -- on its surface where the ion hits).
    -- 0 is a head-on collision; +pi/2 would be a near miss.
    local impact_angle = asin(impact_offset)

    -- In other words, the effect of the above is that impact_angle has
    -- a distribution of p(impact_angle) = sin(2 * impact_angle).

    -- Compute randomized angle [0, 2*pi] for rotation of collision
    -- plane around radial axis.  The is the angle around the
    -- center of the dart board.
    -- Note: all angles are equally likely to hit.
    -- The effect is that impact_theta has a distribution
    -- of p(impact_theta) = 1/(2*pi).
    local impact_theta = 2*pi*rand()

    -- Compute polar coordinates in current velocity reference frame.
    local speed_ion_r, az_ion_r, el_ion_r = rect3d_to_polar3d(vx, vy, vz)

    -- Compute ion velocity components (mm/us).
    local vr_ion = speed_ion_r * cos(impact_angle)    -- radial velocity
    local vt_ion = speed_ion_r * sin(impact_angle)    -- normal velocity

    -- Attenuate ion velocity due to elastic collision.
    -- This is the standard equation for a one-dimensional
    -- elastic collision, assuming the other particle is initially at rest
    -- (in the current reference frame).
    -- Note that the force acts only in the radial direction, which is
    -- normal to the surfaces at the point of contact.
    local vr_ion2 = (vr_ion * (ion_mass - _gas_mass_amu))
                  / (ion_mass + _gas_mass_amu)

    -- Rotate velocity reference frame so that original ion velocity
    -- vector is on the +y axis.
    -- Note: The angle of the new velocity vector with respect to the
    -- +y axis then represents the deflection angle.
    vx, vy, vz = elevation_rotate(90 - deg(impact_angle), vr_ion2, vt_ion, 0)

    -- Rotate velocity reference frame around +y axis.
    -- This rotates the deflection angle and in effect selects the
    -- randomized impact_theta.
    vx, vy, vz = azimuth_rotate(deg(impact_theta), vx, vy, vz)

    -- Rotate velocity reference frame back to the original.
    -- For the incident ion velocity, this would have the effect
    -- of restoring it.
    vx, vy, vz = elevation_rotate(-90 + el_ion_r, vx, vy, vz)
    vx, vy, vz = azimuth_rotate(az_ion_r, vx, vy, vz)

    -- Translate velocity reference frame back to original.
    -- This undoes the prior two translations that make velocity
    -- relative to the colliding gas particle.
    vx = vx + vx_gas + _vx_bar_gas_mmusec
    vy = vy + vy_gas + _vy_bar_gas_mmusec
    vz = vz + vz_gas + _vz_bar_gas_mmusec


    if ion_pz_mm > 180 and ion_pz_mm < 602 then
      -- Set new velocity vector of deflected ion.
      ion_vx_mm, ion_vy_mm, ion_vz_mm = vx, vy, vz
      -- Now lets compute some statistics...

      -- Compute running average of KE.  This is for statistical reporting only.
      -- At thermal equilibrium, KE of the ion and KE of the gas would
      -- be approximately equal according to theory.
      if _trace_level >= 1 then
          -- Compute new ion speed and KE.
          local speed_ion2 = sqrt(ion_vx_mm^2 + ion_vy_mm^2 + ion_vz_mm^2)
          local ke2_ion = speed_to_ke(speed_ion2, ion_mass)

          -- To average ion KE somewhat reliably, we do a running (exponential decay)
          -- average of ion KE over time.  The reset time of the exponential decay
          -- is set to some fraction of the total time-of-flight, so the average
          -- will become more steady as the run progresses (we assume this is a
          -- system that approaches equilibrium).
          -- Note: exp(-x) can be approximated with 1-x for small x.

          -- time between most recent collisions
          local dt = ion_time_of_flight - (last_collision_times[ion_number] or 0)
          -- average over some fraction of TOF
          reset_time = ion_time_of_flight * 0.5
          -- weight for averaging.
          local w = 1 - (dt / reset_time)  -- ~= exp(-dt / reset_time)
          -- update average ion KE
          ke_averages[ion_number] = w * (ke_averages[ion_number] or ke2_ion)
                                  + (1-w) * ke2_ion
          if _trace_level >= 2 then -- more detail
              local T_ion = ke_averages[ion_number] / eV_J / (1.5 * k)
              if trace_count % _trace_skip == 0 then
                  print(string.format(
                      "n=,%d,TOF=,%0.3g,ion KE (eV)=,%0.3e,ion mean KE (eV)=," ..
                      "%0.3e,ion mean temp (K)=,%0.3e",
                      ion_number, ion_time_of_flight, ke2_ion,
                      ke_averages[ion_number], T_ion))
              end
              trace_count = (trace_count + 1) % _trace_skip
          end
          last_collision_times[ion_number] = ion_time_of_flight
      end

      if _mark_collisions ~= 0 then
          mark() -- draw dot at collision point
      end
    end
end

-- SIMION terminate segment. Called on particle termination.
function segment.terminate()
    -- Display some statistics.
    -- Note: At equilibrium, the ion and gas KE become roughly equal.
    if _trace_level >= 1 then
        -- ion temperature
        local T_ion = ke_averages[ion_number] / eV_J / (1.5 * k)
        print(string.format(
            "n=,%d,TOF=,%0.3g,ion mean KE (eV)=,%0.3e,ion mean temp (K)=,%0.3e",
            ion_number, ion_time_of_flight, ke_averages[ion_number], T_ion))
    end
end

function segment.fast_adjust()
 adj_elect06 = ac_voltage * (sin(omega * ion_time_of_flight)) + dc_voltage
 adj_elect07 = ac_voltage * (cos(omega * ion_time_of_flight)) + dc_voltage
 adj_elect08 = ac_voltage * (sin(omega * ion_time_of_flight)) + dc_voltage
 adj_elect09 = ac_voltage * (cos(omega * ion_time_of_flight)) + dc_voltage
 adj_elect10 = ac_voltage * (sin(omega * ion_time_of_flight)) + dc_voltage
 adj_elect11 = ac_voltage * (cos(omega * ion_time_of_flight)) + dc_voltage
end
