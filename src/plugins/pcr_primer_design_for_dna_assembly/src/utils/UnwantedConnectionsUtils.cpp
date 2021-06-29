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
#include <U2Core/PrimerStatistics.h>

namespace U2 {

bool UnwantedConnectionsUtils::hasUnwantedConnections(const QByteArray& forwardSequence,
                                                      double unwantedDeltaG,
                                                      double unwantedMeltingTemperatur,
                                                      int unwantedDimerLength) {
    /**
     * It's reverse complement representation.
     */
    QByteArray revComSequence = DNASequenceUtils::reverseComplement(forwardSequence);

    //Self-dimer
    bool isSelfDimer = isUnwantedSelfDimer(forwardSequence, unwantedDeltaG, unwantedMeltingTemperatur, unwantedDimerLength);

    //TODO:

    return isSelfDimer;
}

bool UnwantedConnectionsUtils::isUnwantedSelfDimer(const QByteArray& forwardSequence,
                                                   double unwantedDeltaG,
                                                   double unwantedMeltingTemperature,
                                                   int unwantedDimerLength) {
    PrimerStatisticsCalculator calc(forwardSequence, PrimerStatisticsCalculator::Direction::DoesntMatter);
    auto dimersInfo = calc.getDimersInfo();
    if (dimersInfo.dimersOverlap.isEmpty()) { // Self dimers aren't found
        return false;
    }

    double dimerMeltingTemp = PrimerStatistics::getMeltingTemperature(dimersInfo.dimer.toLocal8Bit());
    int dimerLength = dimersInfo.dimer.length();
    bool isDeltaGUnwanted = dimersInfo.deltaG < unwantedDeltaG;
    bool isMeltingTemperatureUnwanted = unwantedMeltingTemperature < dimerMeltingTemp;
    bool isLengthUnwanted = unwantedDimerLength < dimerLength;

    return isDeltaGUnwanted && isMeltingTemperatureUnwanted && isLengthUnwanted;
}

bool UnwantedConnectionsUtils::isUnwantedHeteroDimer(const QByteArray& forwardSequence,
                                                     const QByteArray& reverseSequence,
                                                     double unwantedDeltaG,
                                                     double unwantedMeltingTemperatur,
                                                     int unwantedDimerLength) {
    PrimersPairStatistics calc(forwardSequence, reverseSequence);
    auto dimersInfo = calc.getDimersInfo();
    return false;
}


}