__PATCH_ID__ = "DIMJP-2025-12-30-2137-00-tawhiri-models-shape-ascent-profile"

"""
Provide all the balloon models, termination conditions and
functions to combine models and termination conditions.
"""

import calendar
import math
import os

from . import interpolate


_PI_180 = math.pi / 180.0
_180_PI = 180.0 / math.pi


## Up/Down Models #############################################################


def make_constant_ascent(ascent_rate):
    """Return a constant-ascent model at `ascent_rate` (m/s)"""
    def constant_ascent(t, lat, lng, alt):
        return 0.0, 0.0, ascent_rate
    return constant_ascent

def make_piecewise_ascent(ascent_rate, burst_altitude,
                          x1=0.35, f0=1.08, f1=1.31, f2=0.86,
                          min_factor=0.2, steps=2000):
    """
    Piecewise-linear factor:
      (0, f0) -> (x1, f1) -> (1, f2)
    then time-preserving normalization by scale = mean(1/factor)
    """
    def _clip01(x):
        if x < 0.0: return 0.0
        if x > 1.0: return 1.0
        return x

    if burst_altitude <= 0:
        return make_constant_ascent(ascent_rate)

    def raw_factor(x):
        x = _clip01(x)
        if x <= x1:
            # 0..x1
            f = f0 + (f1 - f0) * (x / x1 if x1 > 0 else 0.0)
        else:
            # x1..1
            f = f1 + (f2 - f1) * ((x - x1) / (1.0 - x1) if x1 < 1 else 0.0)
        if f < min_factor:
            f = min_factor
        return f

    # time-preserving scale
    s = 0.0
    for i in range(steps):
        x = (i + 0.5) / steps
        f = raw_factor(x)
        s += 1.0 / f
    scale = s / steps

    def ascent(t, lat, lng, alt):
        x = alt / float(burst_altitude)
        f = raw_factor(x)
        return 0.0, 0.0, ascent_rate * f * scale

    return ascent

def make_piecewise2_ascent(ascent_rate, burst_altitude,
                           x1=0.32, x2=0.75,
                           f0=1.25, f1=1.47, f2=1.05, f3=0.60,
                           min_factor=0.2, steps=4000):
    """
    (0,f0)->(x1,f1)->(x2,f2)->(1,f3)
    time-preserving normalization
    """
    def _clip01(x): return 0.0 if x < 0.0 else (1.0 if x > 1.0 else x)

    if burst_altitude <= 0:
        return make_constant_ascent(ascent_rate)

    def raw_factor(x):
        x = _clip01(x)
        if x <= x1:
            f = f0 + (f1 - f0) * (x / x1 if x1 > 0 else 0.0)
        elif x <= x2:
            f = f1 + (f2 - f1) * ((x - x1) / (x2 - x1) if x2 > x1 else 0.0)
        else:
            f = f2 + (f3 - f2) * ((x - x2) / (1.0 - x2) if x2 < 1 else 0.0)
        return max(f, min_factor)

    s = 0.0
    for i in range(steps):
        x = (i + 0.5) / steps
        s += 1.0 / raw_factor(x)
    scale = s / steps

    def ascent(t, lat, lng, alt):
        x = alt / float(burst_altitude)
        f = raw_factor(x)
        return 0.0, 0.0, ascent_rate * f * scale

    return ascent

def make_alt_dependent_ascent(ascent_rate, burst_altitude,
                              alt_knots, v_knots,
                              min_v=0.2, max_v=20.0,
                              normalize_to_ascent_rate=True,
                              steps=4000):
    """
    Altitude-dependent ascent model: dalt = v(alt)

    alt_knots: 高度[m] の昇順リスト（例: [0, 5000, ..., burst_altitude]）
    v_knots  : その高度での上昇速度[m/s] のリスト（同じ長さ）

    normalize_to_ascent_rate=True:
      v(alt) をスケールして、0→burst の所要時間が
      burst_altitude/ascent_rate になるように補正（=平均上昇率を保つ）
    """

    if burst_altitude <= 0:
        return make_constant_ascent(ascent_rate)

    if len(alt_knots) != len(v_knots) or len(alt_knots) < 2:
        raise ValueError("alt_knots と v_knots は同じ長さ(>=2)にしてください")

    xs = [float(a) for a in alt_knots]
    ys = [float(v) for v in v_knots]

    for i in range(1, len(xs)):
        if xs[i] <= xs[i-1]:
            raise ValueError("alt_knots は昇順にしてください")

    def v_of_alt(alt):
        a = float(alt)

        # 範囲外は端でクランプ
        if a <= xs[0]:
            v = ys[0]
        elif a >= xs[-1]:
            v = ys[-1]
        else:
            # xs[i-1] < a <= xs[i] を探して線形補間
            i = 1
            while xs[i] < a:
                i += 1
            x0, x1 = xs[i-1], xs[i]
            y0, y1 = ys[i-1], ys[i]
            r = (a - x0) / (x1 - x0)
            v = y0 + (y1 - y0) * r

        # 安全クリップ
        if v < min_v:
            v = min_v
        if v > max_v:
            v = max_v
        return v

    # v(alt) を平均上昇率(ascent_rate)に合わせるスケール
    scale = 1.0
    if normalize_to_ascent_rate and ascent_rate > 0:
        # T = ∫ 1/v(alt) dalt を数値積分
        s = 0.0
        for i in range(steps):
            a_mid = (i + 0.5) / steps * float(burst_altitude)
            s += 1.0 / v_of_alt(a_mid)
        T_profile = s / steps * float(burst_altitude)
        T_target = float(burst_altitude) / float(ascent_rate)
        if T_target > 0:
            # v_new = v * scale => T_new = T_profile / scale
            scale = T_profile / T_target

    def ascent(t, lat, lng, alt):
        return 0.0, 0.0, v_of_alt(alt) * scale

    return ascent




def make_drag_descent(sea_level_descent_rate):
    """Return a descent-under-parachute model with sea level descent
       `sea_level_descent_rate` (m/s). Descent rate at altitude is determined
       using an altitude model courtesy of NASA:
       http://www.grc.nasa.gov/WWW/K-12/airplane/atmosmet.html

       For a given altitude the air density is computed, a drag coefficient is
       estimated from the sea level descent rate, and the resulting terminal
       velocity is computed by the returned model function.
    """
    def density(alt):
        temp = pressure = 0.0
        if alt > 25000:
            temp = -131.21 + 0.00299 * alt
            pressure = 2.488 * ((temp + 273.1)/(216.6)) ** (-11.388)
        elif 11000 < alt <= 25000:
            temp = -56.46
            pressure = 22.65 * math.exp(1.73 - 0.000157 * alt)
        else:
            temp = 15.04 - 0.00649 * alt
            pressure = 101.29 * ((temp + 273.1)/288.08) ** (5.256)
        return pressure / (0.2869*(temp + 273.1))

    drag_coefficient = sea_level_descent_rate * 1.1045

    def drag_descent(t, lat, lng, alt):
        return 0.0, 0.0, -drag_coefficient/math.sqrt(density(alt))
    return drag_descent


## Sideways Models ############################################################


def make_wind_velocity(dataset, warningcounts):
    """Return a wind-velocity model, which gives lateral movement at
       the wind velocity for the current time, latitude, longitude and
       altitude. The `dataset` argument is the wind dataset in use.
    """
    get_wind = interpolate.make_interpolator(dataset, warningcounts)
    dataset_epoch = calendar.timegm(dataset.ds_time.timetuple())
    def wind_velocity(t, lat, lng, alt):
        t -= dataset_epoch
        u, v = get_wind(t / 3600.0, lat, lng, alt)
        R = 6371009 + alt
        dlat = _180_PI * v / R
        dlng = _180_PI * u / (R * math.cos(lat * _PI_180))
        return dlat, dlng, 0.0
    return wind_velocity

def _interp1(x, xs, ys):
    """
    1D piecewise-linear interpolation.
    xs: sorted list, same length as ys
    """
    if not xs:
        return 0.0
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    # linear search (xs is short; OK). If long, replace with bisect.
    for i in range(1, len(xs)):
        if x <= xs[i]:
            x0, x1 = xs[i-1], xs[i]
            y0, y1 = ys[i-1], ys[i]
            if x1 == x0:
                return y0
            a = (x - x0) / (x1 - x0)
            return y0 + a * (y1 - y0)
    return ys[-1]


def _clamp(x, lo, hi):
    return lo if x < lo else (hi if x > hi else x)

def _ramp01(alt, a0, a1):
    """
    alt<=a0 -> 0, alt>=a1 -> 1, linear between.
    a0/a1 が None のときは 1 固定（常に適用）
    """
    if a0 is None or a1 is None:
        return 1.0
    if a1 <= a0:
        return 1.0
    if alt <= a0:
        return 0.0
    if alt >= a1:
        return 1.0
    return (alt - a0) / (a1 - a0)
def _clip(x, lo, hi):
    return lo if x < lo else (hi if x > hi else x)

def make_wind_velocity_alt_bias(dataset, warningcounts,
                                alt_knots=None,
                                du_knots=None, dv_knots=None,
                                u_scale_knots=None, v_scale_knots=None,
                                gain_u=1.0, gain_v=1.0,
                                max_delta_u=3.0, max_delta_v=3.0):
    """
    Wind model with altitude-dependent correction, with stabilizers.

    Base wind: (u0,v0) from dataset [m/s]
    Raw corrected:
      u_raw = u0 * su(alt) + du(alt)
      v_raw = v0 * sv(alt) + dv(alt)

    Then apply gains around the base wind:
      u = u0 + gain_u * (u_raw - u0)
      v = v0 + gain_v * (v_raw - v0)

    And clip the deltas:
      (u_raw - u0) is clipped to [-max_delta_u, +max_delta_u]
      (v_raw - v0) is clipped to [-max_delta_v, +max_delta_v]
    """
    alt_knots = alt_knots or []
    du_knots = du_knots or ([0.0] * len(alt_knots))
    dv_knots = dv_knots or ([0.0] * len(alt_knots))
    u_scale_knots = u_scale_knots or ([1.0] * len(alt_knots))
    v_scale_knots = v_scale_knots or ([1.0] * len(alt_knots))

    if not (len(alt_knots) == len(du_knots) == len(dv_knots) ==
            len(u_scale_knots) == len(v_scale_knots)):
        raise ValueError("alt_knots / du_knots / dv_knots / scale_knots length mismatch")

    get_wind = interpolate.make_interpolator(dataset, warningcounts)
    dataset_epoch = calendar.timegm(dataset.ds_time.timetuple())

    def wind_velocity(t, lat, lng, alt):
        t -= dataset_epoch
        u0, v0 = get_wind(t / 3600.0, lat, lng, alt)

        u, v = u0, v0
        if alt_knots:
            su = _interp1(alt, alt_knots, u_scale_knots)
            sv = _interp1(alt, alt_knots, v_scale_knots)
            du = _interp1(alt, alt_knots, du_knots)
            dv = _interp1(alt, alt_knots, dv_knots)

            u_raw = u0 * su + du
            v_raw = v0 * sv + dv

            du_total = u_raw - u0
            dv_total = v_raw - v0

            if max_delta_u is not None:
                du_total = _clip(du_total, -float(max_delta_u), float(max_delta_u))
            if max_delta_v is not None:
                dv_total = _clip(dv_total, -float(max_delta_v), float(max_delta_v))

            u = u0 + float(gain_u) * du_total
            v = v0 + float(gain_v) * dv_total

        R = 6371009 + alt
        dlat = _180_PI * v / R
        dlng = _180_PI * u / (R * math.cos(lat * _PI_180))
        return dlat, dlng, 0.0

    return wind_velocity



def make_reverse_wind_velocity(dataset, warningcounts):
    """Return a reverse wind-velocity model, which gives reverse lateral movement at
       the wind velocity for the current time, latitude, longitude and
       altitude. The `dataset` argument is the wind dataset in use.
       This allows working estimation of a radiosonde's launch site.
    """
    get_wind = interpolate.make_interpolator(dataset, warningcounts)
    dataset_epoch = calendar.timegm(dataset.ds_time.timetuple())
    def wind_velocity(t, lat, lng, alt):
        t -= dataset_epoch
        u, v = get_wind(t / 3600.0, lat, lng, alt)
        # Reverse the sign of the u & v wind components
        u = -1 * u
        v = -1 * v
        R = 6371009 + alt
        dlat = _180_PI * v / R
        dlng = _180_PI * u / (R * math.cos(lat * _PI_180))
        return dlat, dlng, 0.0
    return wind_velocity


## Termination Criteria #######################################################


def make_burst_termination(burst_altitude):
    """Return a burst-termination criteria, which terminates integration
       when the altitude reaches `burst_altitude`.
    """
    def burst_termination(t, lat, lng, alt):
        if alt >= burst_altitude:
            return True
    return burst_termination


def sea_level_termination(t, lat, lng, alt):
    """A termination criteria which terminates integration when
       the altitude is less than (or equal to) zero.

       Note that this is not a model factory.
    """
    if alt <= 0:
        return True

def make_elevation_data_termination(dataset=None):
    """A termination criteria which terminates integration when the
       altitude goes below ground level, using the elevation data
       in `dataset` (which should be a ruaumoko.Dataset).
    """
    def tc(t, lat, lng, alt):
        return (dataset.get(lat, lng) > alt) or (alt <= 0)
    return tc

def make_time_termination(max_time):
    """A time based termination criteria, which terminates integration when
       the current time is greater than `max_time` (a UNIX timestamp).
    """
    def time_termination(t, lat, lng, alt):
        if t > max_time:
            return True
    return time_termination

def make_dummy_termination():
    """A dummy termination criteria, which immediately terminates """
    def dum(t, lat, lng, alt):
        return True
    return dum


## Model Combinations #########################################################


def make_linear_model(models):
    """Return a model that returns the sum of all the models in `models`.
    """
    def linear_model(t, lat, lng, alt):
        dlat, dlng, dalt = 0.0, 0.0, 0.0
        for model in models:
            d = model(t, lat, lng, alt)
            dlat, dlng, dalt = dlat + d[0], dlng + d[1], dalt + d[2]
        return dlat, dlng, dalt
    return linear_model


def make_any_terminator(terminators):
    """Return a terminator that terminates when any of `terminators` would
       terminate.
    """
    def terminator(t, lat, lng, alt):
        return any(term(t, lat, lng, alt) for term in terminators)
    return terminator


## Pre-Defined Profiles #######################################################


def standard_profile(ascent_rate, burst_altitude, descent_rate,
                     wind_dataset, elevation_dataset, warningcounts):

    # === paste from fit_wind_bias_knots.py ===
    ALT_KNOTS = [1000, 3000, 5000, 7000, 9000, 11000, 13000, 15000, 17000, 19000, 21000, 23000, 25000]
    DU_KNOTS  = [1.140, 1.795, 1.795, 0.357, 0.357, 0.674, 0.674, -3.870, -3.870, -1.384, -1.094, -1.384, -1.555]  # +east m/s
    DV_KNOTS  = [-0.347, -0.719, -0.719, -1.200, -1.200, -2.550, 3.311, 3.311, 3.311, -0.240, -0.240, 0.732, 1.556]  # +north m/s
    # ===========================

    SU_KNOTS = [1.0] * len(ALT_KNOTS)
    SV_KNOTS = [1.0] * len(ALT_KNOTS)

    # === stabilizers / gains ===
    # まず弱めに開始（ここをいじってチューニングする）
    gain_u_up   = 0.4
    gain_v_up   = 0.4
    gain_u_down = 0.4
    gain_v_down = 0.0   # ★まずは降下の緯度悪化を止めるため 0 推奨

    max_du = 2.0        # ★最初は 1〜2 m/s 程度で
    max_dv = 2.0

    wind_up = make_wind_velocity_alt_bias(
        wind_dataset, warningcounts,
        alt_knots=ALT_KNOTS,
        du_knots=DU_KNOTS, dv_knots=DV_KNOTS,
        u_scale_knots=SU_KNOTS, v_scale_knots=SV_KNOTS,
        gain_u=gain_u_up, gain_v=gain_v_up,
        max_delta_u=max_du, max_delta_v=max_dv
    )

    wind_down = make_wind_velocity_alt_bias(
        wind_dataset, warningcounts,
        alt_knots=ALT_KNOTS,
        du_knots=DU_KNOTS, dv_knots=DV_KNOTS,
        u_scale_knots=SU_KNOTS, v_scale_knots=SV_KNOTS,
        gain_u=gain_u_down, gain_v=gain_v_down,  # ★v補正を切れる
        max_delta_u=max_du, max_delta_v=max_dv
    )

    model_up = make_linear_model([
        make_constant_ascent(ascent_rate),
        wind_up
    ])
    term_up = make_burst_termination(burst_altitude)

    model_down = make_linear_model([
        make_drag_descent(descent_rate),
        wind_down
    ])
    term_down = make_elevation_data_termination(elevation_dataset)

    return ((model_up, term_up), (model_down, term_down))




def float_profile(ascent_rate, float_altitude, stop_time, dataset, warningcounts):
    """Make a model chain for the typical floating balloon situation of ascent
       at constant altitude to a float altitude which persists for some
       amount of time before stopping. Descent is in general not modelled.
    """

    model_up = make_linear_model([make_constant_ascent(ascent_rate),
                                  make_wind_velocity(dataset, warningcounts)])
    term_up = make_burst_termination(float_altitude)
    model_float = make_wind_velocity(dataset, warningcounts)
    term_float = make_time_termination(stop_time)

    return ((model_up, term_up), (model_float, term_float))


def reverse_profile(ascent_rate, wind_dataset, elevation_dataset, warningcounts):
    """Make a model chain used to estimate a balloon's launch site location, based on
       the current position, and a known ascent rate. This model only works for a balloon
       on ascent.

       Requires the balloon `ascent_rate`,
       and additionally requires the dataset to use for wind velocities.

       Returns a tuple of (model, terminator) pairs.
    """

    model_up = make_linear_model([make_constant_ascent(ascent_rate),
                                  make_wind_velocity(wind_dataset, warningcounts)])

    term_up = make_dummy_termination()

    model_down = make_linear_model([make_constant_ascent(abs(ascent_rate)),
                                    make_wind_velocity(wind_dataset, warningcounts)])
    term_down = make_elevation_data_termination(elevation_dataset)

    return ((model_up, term_up), (model_down, term_down))