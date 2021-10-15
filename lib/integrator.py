#!/usr/bin/python
# -*- coding:utf-8 -*-
"""
Implementation of the various integrators for numerical integration.

Comes from the assumption that the problem is analytically defined in position-momentum (q-p) space for a given hamiltonian H.
"""
import numpy as np

def dp_dt(m_array, q_array):
    """
    Time derivative of the momentum, given by the position derivative of the Hamiltonian.
    dp/dt = -dH/dq
    """
    dp_array = np.zeros(q_array.shape)
    for i in range(q_array.shape[0]):
        m_j = np.delete(m_array, i)
        q_j = np.delete(q_array, i, 0)
        dp_array = m_array[i]*np.sum((m_j/np.sum((q_j-q_array[i])**3, axis=1)).reshape((q_j.shape[0],1))*(q_j-q_array[i]), axis=0)
    return dp_array

def frogleap(duration, step, m_array, q_array, p_array):
    """
    Leapfrog integrator for first order partial differential equations.
    iteration : half-step drift -> full-step kick -> half-step drift
    """
    N = np.ceil(duration/step).astype(int)
    for _ in range(N):
        # half-step drift
        q_array, p_array = q_array + step/2*p_array/m_array , p_array
        # full-step kick
        q_array, p_array = q_array , p_array - step*dp_dt(m_array, q_array)
        # half-step drift
        q_array, p_array = q_array + step/2*p_array/m_array , p_array

    return q_array, p_array
