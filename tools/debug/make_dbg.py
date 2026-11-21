# make_dbg.py -- diagnostic copy of the MPI solver: with -nlogit 2, every Newton iteration also prints one 'DBG' line with the
# number of elements on each non-smooth switch of the local update (drag bound per slip system, GB-spin cap, slip-direction sign
# changes since the previous iteration, the merit accepted by the line search and the merit after the repeated local update, elements whose cap/bound state changed, dual-flux steps cut at |w| = 1, elements with |w| ~ 1) and the state of the element with the largest |abar|.
#   python3 make_dbg.py [cosserat_kwc_v2_mpi.edp] [cosserat_kwc_v2_dbg_mpi.edp]
import sys
src = sys.argv[1] if len(sys.argv) > 1 else 'cosserat_kwc_v2_mpi.edp'
dst = sys.argv[2] if len(sys.argv) > 2 else 'cosserat_kwc_v2_dbg_mpi.edp'
s = open(src).read()
def rep(a, b):
    global s
    n = s.count(a); assert n == 1, ('target not found once', a[:80], n); s = s.replace(a, b)
a = 'int kwcTan = atoi(arg("-kwctan", "0")), almC = atoi(arg("-almc", "0")), nlogIt = atoi(arg("-nlogit", "0"));'
rep(a, a + '\nreal[int] dbgSg(3*nt), dbgSt(4*nt); dbgSg = 0; dbgSt = 0; real dbgSwN = 0;   // make_dbg.py: previous slip signs and switch states; truncated dual steps')
a = '          wdx[k] += sw*ddx[k]; wdy[k] += sw*ddy[k];'
rep(a, a + ' if (sw < 1 && ownE[k] > 0) dbgSwN += 1;   // make_dbg.py')
a = '    if (nlogIt) nlog << curStep << "," << it << "," << rn << "," << clock() << "," << aMax << "," << betaIt*180/pi << endl;'
rep(a, a + r'''
    if (nlogIt > 1) {   // make_dbg.py
      real[int] st(4*nt), sg(3*nt), fl(nt), dd(nt);
      dragBound(ta, eEa); fl = abs(dpa); dd = dpmax; dd *= (1 - 1e-9); fl -= dd; for [i, v : fl] v = (v >= 0); st(0:nt-1) = fl;
      dragBound(tb, eEb); fl = abs(dpb); dd = dpmax; dd *= (1 - 1e-9); fl -= dd; for [i, v : fl] v = (v >= 0); st(nt:2*nt-1) = fl;
      dragBound(tcc, eEc); fl = abs(dpc); dd = dpmax; dd *= (1 - 1e-9); fl -= dd; for [i, v : fl] v = (v >= 0); st(2*nt:3*nt-1) = fl;
      fl = abs(dwst); fl -= spinCap*(1 - 1e-9); for [i, v : fl] v = (spinCap > 0 && v >= 0); st(3*nt:4*nt-1) = fl;
      for [i, v : fl] v = (ta[i] > 0) - (ta[i] < 0); sg(0:nt-1) = fl;
      for [i, v : fl] v = (tb[i] > 0) - (tb[i] < 0); sg(nt:2*nt-1) = fl;
      for [i, v : fl] v = (tcc[i] > 0) - (tcc[i] < 0); sg(2*nt:3*nt-1) = fl;
      real[int] cnt(6); cnt = 0;
      for (int e = 0; e < nt; e++) if (ownE[e] > 0) {
        cnt[0] += st[e]; cnt[1] += st[nt+e]; cnt[2] += st[2*nt+e]; cnt[3] += st[3*nt+e];
        if (it > 0) { bool f = 0, c = 0; for (int k = 0; k < 3; k++) f = f || (sg[k*nt+e] != dbgSg[k*nt+e]); for (int k = 0; k < 4; k++) c = c || (st[k*nt+e] != dbgSt[k*nt+e]); cnt[4] += f; cnt[5] += c; }
      }
      real[int] cg(6); mpiAllReduce(cnt, cg, mpiCommWorld, mpiSUM);
      real wn1 = 0, swg = 0; { real c1 = 0; for (int e = 0; e < nt; e++) if (ownE[e] > 0 && wdx[e]^2 + wdy[e]^2 > 0.9999^2) c1 += 1; mpiAllReduce(c1, wn1, mpiCommWorld, mpiSUM); mpiAllReduce(dbgSwN, swg, mpiCommWorld, mpiSUM); dbgSwN = 0; }
      dbgSg = sg; dbgSt = st;
      Ph xcE = x, ycE = y;
      real[int] qa = abs(aMean); qa = qa .* ownE; int ie = qa.imax;
      if (mpirank == 0) cout << "DBG " << curStep << " it " << it << ": drag-bound elements a/b/c " << cg[0] << " " << cg[1] << " " << cg[2] << ", spin-capped " << cg[3]
                             << ", slip-sign changes " << cg[4] << ", elements with a changed bound/cap state " << cg[5] << ", dual steps truncated at |w| = 1 (previous update) " << swg << ", |w| > 0.9999 " << wn1 << endl;
      if (qa[ie] >= aMax*(1 - 1e-12) && aMax > 0) cout << "DBG " << curStep << " it " << it << " max|abar| element " << ie << " (rank " << mpirank << ") at (" << xcE[][ie] << ", " << ycE[][ie] << "): abar " << aMean[ie]
           << " lambda " << lamC[ie] << " dgamma " << dpa[ie] << " " << dpb[ie] << " " << dpc[ie] << " dw* " << dwst[ie] << " tau " << ta[ie] << " " << tb[ie] << " " << tcc[ie]
           << " tauc " << tauc[ie] << " state " << st[ie] << st[nt+ie] << st[2*nt+ie] << st[3*nt+ie] << endl;
    }''')
a = '        if (rt == rt && rt <= (1 - 1e-4*w)*rn) { lsOk = 1; break; }'
rep(a, '        if (rt == rt && rt <= (1 - 1e-4*w)*rn) { lsOk = 1; if (nlogIt > 1 && mpirank == 0) cout << "DBG " << curStep << " it " << it << " line search accepted w " << w << " merit " << rt << endl; break; }   // make_dbg.py')
a = '    mark(tg + "local plastic update (slip/GB-spin Newton per element)");\n    localUpdate();\n'
rep(a, a + '    if (nlogIt > 1 && it > 0) { residual(); real r1 = sqrt(gdotW(rhs, rhs)); if (almC == 4) r1 = sqrt(r1*r1 + conRes()); real[int] q1 = abs(aMean); real a1 = gmaxS(q1.max); if (mpirank == 0) cout << "DBG " << curStep << " it " << it << " after the repeated local update (before the tangent): merit " << r1 << " max|abar| " << a1 << endl;\n      { real[int] rq(nt); localRes(rq); int iq = rq.imax; real lm = rq[iq], lg = gmaxS(lm); real nb = 0; for (int e = 0; e < nt; e++) if (ownE[e] > 0 && rq[e] > 1e-10) nb += 1; real nbg = nb; mpiAllReduce(nb, nbg, mpiCommWorld, mpiSUM);\n        if (mpirank == 0) cout << "DBG " << curStep << " it " << it << " local residual after the update: max " << lg << ", elements above 1e-10: " << nbg << endl;\n        if (lm >= lg*(1 - 1e-12) && lg > 1e-10) cout << "DBG " << curStep << " it " << it << " worst local element " << iq << " (rank " << mpirank << "): r " << ra[iq] << " " << rb[iq] << " " << rc[iq] << " " << rwl[iq] << " dgamma " << dpa[iq] << " " << dpb[iq] << " " << dpc[iq] << " dw* " << dwst[iq] << " tau " << ta[iq] << " " << tb[iq] << " " << tcc[iq] << " tau* " << tstar[iq] << " tauc " << tauc[iq] << " phi-weight kspin " << kspin[iq] << endl; residual(); }   // make_dbg.py\n    }\n')
open(dst, 'w').write(s); print('wrote', dst)
