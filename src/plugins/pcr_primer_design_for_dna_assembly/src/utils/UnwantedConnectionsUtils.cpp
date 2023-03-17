/**
 * UGENE - Integrated Bioinformatics Tools.
 * Copyright (C) 2008-2021 UniPro <ugene@unipro.ru>
 * http://ugene.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */

#include "UnwantedConnectionsUtils.h"

#include <U2Core/DNASequenceUtils.h>
#include <U2Core/Log.h>
#include <U2Core/PrimerStatistics.h>

#include "PCRPrimerDesignForDNAAssemblyPlugin.h"

namespace U2 {

bool UnwantedConnectionsUtils::isUnwantedSelfDimer(const QByteArray& forwardSequence,
                                                   double minGibbs,
                                                   double maxTemp,
                                                   int maxLenth,
                                                   const QSharedPointer<TmCalculator>& tmCalculator) {
    QString unused;
    return isUnwantedSelfDimer(forwardSequence, minGibbs, maxTemp, maxLenth, tmCalculator,
                               unused);
}

bool UnwantedConnectionsUtils::isUnwantedSelfDimer(const QByteArray &forwardSequence,
                                                   double minGibbs,
                                                   double maxTemp,
                                                   int maxLenth,
                                                   const QSharedPointer<TmCalculator>& tmCalculator,
                                                   QString &report) {
    PrimerStatisticsCalculator calc(forwardSequence, tmCalculator, PrimerStatisticsCalculator::Direction::DoesntMatter);
    report = PCRPrimerDesignForDNAAssemblyPlugin::tr("<b>Self-dimer:</b><br>");
    return areUnwantedParametersPresentedInDimersInfo(calc.getDimersInfo(), minGibbs, maxTemp,
                                                      maxLenth, tmCalculator, report);
}

bool UnwantedConnectionsUtils::isUnwantedHeteroDimer(const QByteArray &forwardSequence,
                                                     const QByteArray &reverseSequence,
                                                     double minGibbs,
                                                     double maxTemp,
                                                     int maxLenth,
                                                     const QSharedPointer<TmCalculator>& tmCalculator) {
    QString unused;
    return isUnwantedHeteroDimer(forwardSequence, reverseSequence, minGibbs, maxTemp,
                                 maxLenth, tmCalculator, unused);
}

bool UnwantedConnectionsUtils::isUnwantedHeteroDimer(const QByteArray& forwardSequence,
                                                     const QByteArray& reverseSequence,
                                                     double minGibbs,
                                                     double maxTemp,
                                                     int maxLenth,
                                                     const QSharedPointer<TmCalculator>& tmCalculator,
                                                     QString &report) {
    PrimersPairStatistics calc(forwardSequence, reverseSequence, tmCalculator);
    report = PCRPrimerDesignForDNAAssemblyPlugin::tr("<b>Hetero-dimer:</b><br>");
    return areUnwantedParametersPresentedInDimersInfo(calc.getDimersInfo(), minGibbs, maxTemp,
                                                      maxLenth, tmCalculator, report);
}

bool UnwantedConnectionsUtils::areUnwantedParametersPresentedInDimersInfo(const DimerFinderResult& dimersInfo,
                                                                          double minGibbs,
                                                                          double maxTemp,
                                                                          int maxLenth,
                                                                          const QSharedPointer<TmCalculator>& tmCalculator,
                                                                          QString &report) {
    if (dimersInfo.dimersOverlap.isEmpty()) {
        return false;
    }
    double dimerMeltingTemp = tmCalculator->getMeltingTemperature(dimersInfo.dimer.toLocal8Bit());
    int dimerLength = dimersInfo.dimer.length();
    bool isDeltaGUnwanted = dimersInfo.deltaG <= minGibbs;
    bool isMeltingTemperatureUnwanted = maxTemp <= dimerMeltingTemp;
    bool isLengthUnwanted = maxLenth <= dimerLength;
    bool isUnwantedParameter = isDeltaGUnwanted && isMeltingTemperatureUnwanted && isLengthUnwanted;
    report += PCRPrimerDesignForDNAAssemblyPlugin::tr(dimersInfo.getFullReport().toLocal8Bit());
    if (isUnwantedParameter) {
        algoLog.details(report);
    }

    return isUnwantedParameter;
}

}  // namespace U2
