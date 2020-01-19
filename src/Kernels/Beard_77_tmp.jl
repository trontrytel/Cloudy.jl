module Beard_77_tmp

export Beard77_sea_level_velocity
export Beard77_table3

function Beard77_table3()

  # drop dimeter [m]
  diameter = Real[1, 1.5, 2, 3, 4, 6, 8, 10, 15, 20, 30, 40, 60, 80, 100, 150,
                  200, 300, 400, 600, 800, 1000, 1500, 2000, 3000, 4000, 5000,
                  6000] * 1e-6

  # terminal velocity [m/s]
  vel_sea = Real[3.48e-3, 7.47e-3, 1.29e-2, 2.84e-2, 4.98e-2, 1.105e-1,
                 1.951e-1, 0.304, 0.670, 1.201, 2.68, 4.70, 10.18, 17.08, 25,
                 46.7, 69.4, 114.7, 158.8, 244, 326, 402, 542, 653, 808, 885,
                 912, 914] * 1e-2

  return [diameter, vel_sea]
end

function Beard77_sea_level_velocity(diameter::FT) where {FT<:Real}

  c_sml = FT[0.105035e2, 0.108750e1, -0.133245, -0.659969e-2]
  c_big = FT[0.65639e1, -0.10391e1, -0.14001e1, -0.82736e0, -0.34277e0,
             -0.83072e-1, -0.10583e-1, -0.54208e-3]

  vel_sea = 44

  if diameter < 1e-6
    vel_sea = FT(0)
  end

  if diameter >= 1e-6 && diameter < 40e-6
    x = log(diameter)
    y = FT(0)
    for i in 0:3
      y += c_sml[i+1] * x^i
    end
    vel_sea = exp(y)
  end

  if diameter >= 40e-6 && diameter < 6e-3
    x = log(diameter)
    y = FT(0)
    for i in 0:7
      y += c_big[i+1] * x^i
    end
    vel_sea = exp(y)
  end

  if diameter >= 6e-3
    x=log(6e-3)
    y = FT(0)
    for i in 0:7
      y += c_big[i+1] * x^i
    end
    vel_sea = exp(y)
  end

  return vel_sea
end

end
